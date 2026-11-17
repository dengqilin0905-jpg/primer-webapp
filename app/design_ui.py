from sequencing_ui import sequencing_download
import json
import re
import time
from typing import Any

import streamlit as st

from db import get_conn
from backbone_presets import BackbonePresetError, list_backbone_presets
from design_config import config_provenance, custom_backbone_design_config, effective_remove_stop, shared_preset_design_config
from excel_export import export_from_template_bytes, validate_primer_rows
from fragment_export import build_fragment_length_simple_txt, build_fragment_length_txt, build_non_cds_amplicon_fasta
from genbank_export import build_genbank_bundle
from primer_design import design_primers_for_gene
from snapgene_export import build_snapgene_dna_bundle, find_snapgene_cli
from templates_ui import get_default_template_id, get_effective_template_bytes, get_user_templates
from template_resolver import GLOBAL_DEFAULT_TEMPLATE_PATH, resolve_order_template
from naming import build_construct_name, build_primer_label, build_map_filename


def _primer_name(backbone_name: str, gene_name: str, direction: str) -> str:
    return build_primer_label(gene_name, direction, backbone=backbone_name)


def _remove_stop_label(seq_type: str, remove_stop: bool) -> str:
    if (seq_type or "").upper() != "CDS":
        return "N/A"
    return "Yes" if remove_stop else "No"


def _bb_label(backbone_row) -> str:
    enz = (backbone_row["enzyme_sites"] or "").strip()
    return f"{backbone_row['name']} [{enz}]" if enz else f"{backbone_row['name']} [no-enzyme]"


def _source_ref(config: dict[str, Any]) -> dict[str, Any]:
    return {
        "backbone_source": config["backbone_source"],
        "backbone_id": config.get("backbone_id"),
        "backbone_name": config["backbone_name"],
        "backbone_label": config["backbone_label"],
        "preset_id": config.get("preset_id"),
        "preset_version": config.get("preset_version"),
        "vector_sha256": config.get("vector_sha256"),
    }


def _resolve_source_config(ref: dict[str, Any], backbones_by_id: dict[int, Any]) -> dict[str, Any]:
    if ref.get("backbone_source") == "shared_preset":
        return shared_preset_design_config(ref["preset_id"], version=ref.get("preset_version"))
    row = backbones_by_id.get(int(ref["backbone_id"]))
    if row is None:
        raise ValueError(f"Backbone no longer exists: {ref.get('backbone_label', ref.get('backbone_id'))}")
    return custom_backbone_design_config(row)


def _normalize_mapping_cart(cart, backbones, source_refs=None):
    by_id = {int(b["id"]): b for b in backbones}
    by_name_first = {}
    for b in backbones:
        by_name_first.setdefault(b["name"], b)

    fixed = []
    for item in cart:
        if item.get("backbone_source") == "shared_preset":
            if not item.get("preset_id"):
                continue
            fixed.append(dict(item))
            continue
        seq_type = str(item.get("seq_type") or "CDS")
        remove_stop = bool(item.get("remove_stop")) if seq_type.upper() == "CDS" else False
        gene_names = list(item.get("gene_names") or [])
        if not gene_names:
            continue

        bb = None
        if "backbone_id" in item and item.get("backbone_id") is not None:
            bb = by_id.get(int(item["backbone_id"]))
        if bb is None and item.get("backbone_name"):
            bb = by_name_first.get(item["backbone_name"])
        if bb is None:
            continue

        fixed.append(
            {
                "backbone_id": int(bb["id"]),
                "backbone_name": bb["name"],
                "backbone_label": _bb_label(bb),
                "gene_names": gene_names,
                "seq_type": seq_type,
                "remove_stop": remove_stop,
                "backbone_source": "user_backbone",
                "preset_id": None,
                "preset_version": None,
                "vector_sha256": None,
            }
        )
    return fixed


def _build_preview_rows(mapping_cart):
    out = []
    for item in mapping_cart:
        for gene in item["gene_names"]:
            out.append(
                {
                    "backbone": item["backbone_label"],
                    "gene": gene,
                    "seq_type": item["seq_type"],
                    "remove_stop": _remove_stop_label(item["seq_type"], bool(item.get("remove_stop"))),
                    "expected_F": _primer_name(item["backbone_name"], gene, "F"),
                    "expected_R": _primer_name(item["backbone_name"], gene, "R"),
                }
            )
    return out


def _flatten_combos(mapping_cart) -> list[dict[str, Any]]:
    combos = []
    for item in mapping_cart:
        seq_type = str(item.get("seq_type") or "CDS")
        remove_stop = bool(item.get("remove_stop")) if seq_type.upper() == "CDS" else False
        for gene in item.get("gene_names", []):
            combos.append(
                {
                    "backbone_id": int(item["backbone_id"]) if item.get("backbone_id") is not None else None,
                    "backbone_name": item["backbone_name"],
                    "backbone_label": item["backbone_label"],
                    "gene_name": gene,
                    "seq_type": seq_type,
                    "remove_stop": remove_stop,
                    "backbone_source": item.get("backbone_source", "user_backbone"),
                    "preset_id": item.get("preset_id"),
                    "preset_version": item.get("preset_version"),
                    "vector_sha256": item.get("vector_sha256"),
                }
            )
    return combos


def _cartesian_items(source_refs, genes, seq_type, remove_stop):
    items = []
    for ref in source_refs:
        items.append(
            {
                "backbone_id": ref.get("backbone_id"),
                "backbone_name": ref["backbone_name"],
                "backbone_label": ref["backbone_label"],
                "gene_names": list(genes),
                "seq_type": seq_type,
                "remove_stop": bool(remove_stop) if seq_type == "CDS" else False,
                "backbone_source": ref["backbone_source"],
                "preset_id": ref.get("preset_id"),
                "preset_version": ref.get("preset_version"),
                "vector_sha256": ref.get("vector_sha256"),
            }
        )
    return items


def _format_batch_mode(x: str) -> str:
    return {"mapping": "Mapping (existing)", "cartesian": "Cartesian (Gene × Backbone)"}.get(x, x)


def _with_map_options(
    rows: list[dict[str, Any]],
    feature_type: str,
    feature_name: str,
    feature_note: str,
    homology_type: str,
) -> list[dict[str, Any]]:
    out = []
    for row in rows:
        item = dict(row)
        item["insert_feature_type"] = feature_type
        item["insert_feature_name"] = (feature_name or "").strip() or item.get("gene_name") or "insert"
        item["insert_feature_note"] = (feature_note or "").strip() or "Inserted sequence"
        item["homology_feature_type"] = homology_type
        out.append(item)
    return out


def design_page(user, config):
    st.subheader("Design (Batch)")
    conn = get_conn()

    prof = conn.execute("SELECT * FROM user_profiles WHERE user_id=?", (user["id"],)).fetchone()
    profile = dict(prof) if prof else None

    backbones = conn.execute(
        "SELECT * FROM backbones WHERE user_id=? ORDER BY name, enzyme_sites, id",
        (user["id"],),
    ).fetchall()
    try:
        shared_presets = list_backbone_presets()
    except BackbonePresetError as exc:
        shared_presets = []
        st.error(f"Shared backbone preset registry unavailable: {exc}")
    shared_refs = []
    for p in shared_presets:
        try:
            shared_refs.append(_source_ref(shared_preset_design_config(p["preset_id"], version=p["preset_version"])))
        except (BackbonePresetError, ValueError) as exc:
            st.error(f"Shared preset {p.get('preset_id', 'unknown')} is unavailable: {exc}")
    custom_refs = [
        _source_ref(custom_backbone_design_config(row))
        for row in backbones
    ]
    source_refs = shared_refs + custom_refs
    if not source_refs:
        st.warning("Please add a custom backbone or use a shared preset.")
        conn.close()
        return

    st.markdown("#### Backbone source")
    source_kind = st.radio(
        "Backbone source",
        (["Shared preset"] if shared_refs else []) + (["My/custom backbone"] if custom_refs else []),
        horizontal=True,
        key="design_backbone_source",
    )
    if source_kind == "Shared preset":
        source_labels = [r["backbone_label"] for r in shared_refs]
        picked_source_label = st.selectbox("Select shared preset", source_labels, key="design_shared_preset_pick")
        selected_source_ref = next(r for r in shared_refs if r["backbone_label"] == picked_source_label)
    else:
        source_labels = [r["backbone_label"] for r in custom_refs]
        picked_source_label = st.selectbox("Select custom backbone", source_labels, key="design_custom_backbone_pick")
        selected_source_ref = next(r for r in custom_refs if r["backbone_label"] == picked_source_label)
    try:
        selected_config = _resolve_source_config(selected_source_ref, {int(b["id"]): b for b in backbones})
    except (BackbonePresetError, ValueError) as exc:
        st.error(f"Selected backbone configuration is invalid: {exc}")
        conn.close()
        return
    if selected_config["backbone_source"] == "shared_preset":
        st.info(
            f"Shared preset {selected_config['preset_id']} v{selected_config['preset_version']} "
            f"| status={selected_config['preset_status']} | validation={selected_config['validation_level']} "
            f"| user-approved={selected_config['approved_by_user']} | {selected_config['enzyme_sites']} "
            f"| homology {selected_config['homology_arm_length']} bp "
            f"| {selected_config['fusion_order']} | remove_stop default={selected_config['remove_stop_default']}"
        )
        st.caption("Fusion-frame validation is required by this preset but is not automatically executed by the primer engine.")
    else:
        st.caption(f"Using selected custom backbone: {selected_config['backbone_label']}")

    defaults = dict(selected_config.get("design_defaults") or {})
    defaults.setdefault("remove_stop", selected_config.get("remove_stop_default", True))

    st.markdown("#### Basic parameters")
    st.caption("Key choices that affect the final primers stay visible here. Advanced chemistry and warning settings are folded below.")
    col1, col2 = st.columns(2)
    with col1:
        default_seq_type = st.selectbox("Default sequence type", ["CDS", "Promoter", "Other"], key="design_default_seq_type")
        default_remove_stop = st.checkbox(
            "Default remove terminal stop codon (CDS only)",
            value=bool(defaults.get("remove_stop", True)),
            key="design_default_remove_stop",
        )
        min_len = st.number_input("Min primer length", value=int(defaults.get("min_len", 18)), step=1, key="design_min_len")
        max_len = st.number_input("Max primer length", value=int(defaults.get("max_len", 30)), step=1, key="design_max_len")
    with col2:
        tm_min = st.number_input("Tm min (degC)", value=float(defaults.get("tm_min", 55.0)), step=0.5, key="design_tm_min")
        tm_max = st.number_input("Tm max (degC)", value=float(defaults.get("tm_max", 65.0)), step=0.5, key="design_tm_max")
        target_tm = st.number_input("Target Tm (degC)", value=float(defaults.get("target_tm", 60.0)), step=0.5, key="design_target_tm")
        tm_pair_delta_max = st.number_input("Max |TmF-TmR| (soft)", value=2.0, min_value=0.0, step=0.5, key="design_tm_pair_delta_max")

    if selected_config["backbone_source"] == "shared_preset" and bool(default_remove_stop) != bool(selected_config["remove_stop_default"]):
        st.caption("Current remove_stop value is an explicit construct override of the selected preset default.")

    with st.expander("Advanced parameters", expanded=False):
        adv1, adv2, adv3 = st.columns(3)
        with adv1:
            gc_min = st.number_input("GC min (%)", value=40.0, min_value=0.0, max_value=100.0, step=1.0, key="design_gc_min")
            gc_max = st.number_input("GC max (%)", value=60.0, min_value=0.0, max_value=100.0, step=1.0, key="design_gc_max")
            non_cds_forbid_homopolymer = st.checkbox("Forbid AAAAA/TTTTT in non-CDS primers", value=True, key="design_non_cds_forbid_hp")
        with adv2:
            tm_priority_weight = st.number_input("Tm priority weight", value=2.5, min_value=0.1, step=0.1, key="design_tm_priority_weight")
            cds_recovery_enable = st.checkbox("Enable CDS rescue when strict filters fail", value=True, key="design_cds_recovery_enable")
            tm_delta_mode = st.selectbox("Tm delta mode", ["soft", "hard"], key="design_tm_delta_mode")
        with adv3:
            enable_frame_check = st.checkbox(
                "Enable primer-engine frame check",
                value=bool(selected_config.get("primer_engine_frame_check", False)),
                disabled=(selected_config["backbone_source"] == "shared_preset"),
                key="design_enable_frame_check",
            )
            construct_type = st.selectbox(
                "Construct type",
                ["CDS_no_tag", "N_terminal_tag", "C_terminal_tag", "promoter_or_non_coding"],
                index=(
                    ["CDS_no_tag", "N_terminal_tag", "C_terminal_tag", "promoter_or_non_coding"].index(selected_config["construct_type"])
                    if selected_config.get("construct_type") in {"CDS_no_tag", "N_terminal_tag", "C_terminal_tag", "promoter_or_non_coding"}
                    else 0
                ),
                disabled=(selected_config["backbone_source"] == "shared_preset"),
                key="design_construct_type",
            )
        tm_profile_label = st.selectbox("Tm Profile", ["SnapGene Preset", "Custom"], index=0, key="design_tm_profile")
        if tm_profile_label == "Custom":
            Na = st.number_input("Na+ (mM)", value=50.0, step=5.0, key="design_na")
            Mg = st.number_input("Mg2+ (mM)", value=1.5, step=0.5, key="design_mg")
            dNTPs = st.number_input("dNTPs (mM)", value=0.2, step=0.1, key="design_dntps")
            primer_conc = st.number_input("Primer conc (nM)", value=250.0, step=50.0, key="design_primer_conc")
            tm_profile = "custom"
        else:
            tm_profile = "snapgene_preset"
            Na, Mg, dNTPs, primer_conc = 50.0, 1.5, 0.2, 250.0
            st.caption("Preset params: Na=50.0 mM, Mg=1.5 mM, dNTPs=0.2 mM, primer=250 nM")

    policy_defaults = {
        "gc_min": float(gc_min),
        "gc_max": float(gc_max),
        "length_mode": "adaptive_soft",
        "cds_len_pref": [18, 21],
        "non_cds_len_pref": [20, 26],
        "tm_priority_weight": float(tm_priority_weight),
        "non_cds_forbid_homopolymer": bool(non_cds_forbid_homopolymer),
        "tm_pair_delta_max": float(tm_pair_delta_max),
        "tm_delta_mode": str(tm_delta_mode),
        "tm_profile": tm_profile,
        "cds_recovery_enable": bool(cds_recovery_enable),
    }
    st.caption(f"Design policy: Tm-priority + GC {gc_min:.0f}-{gc_max:.0f}% + adaptive soft length + {tm_delta_mode} |TmF-TmR| + CDS rescue.")

    search = st.text_input("Filter genes by name (optional)", key="design_gene_filter")
    if search:
        genes = conn.execute(
            "SELECT id, gene_name, sequence FROM genes WHERE user_id=? AND gene_name LIKE ? ORDER BY gene_name LIMIT 5000",
            (user["id"], f"%{search}%"),
        ).fetchall()
    else:
        genes = conn.execute(
            "SELECT id, gene_name, sequence FROM genes WHERE user_id=? ORDER BY gene_name LIMIT 5000",
            (user["id"],),
        ).fetchall()

    batch_mode = st.radio("Batch mode", ["mapping", "cartesian"], format_func=_format_batch_mode, horizontal=True, key="design_batch_mode")

    if "mapping_cart" not in st.session_state:
        st.session_state["mapping_cart"] = []
    st.session_state["mapping_cart"] = _normalize_mapping_cart(st.session_state["mapping_cart"], backbones)

    cart_items = []
    if batch_mode == "mapping":
        st.markdown("#### Build mapping (Backbone -> Genes)")
        picked_label = selected_source_ref["backbone_label"]
        st.caption(f"Selected backbone source: {picked_label}")

        selected_genes = st.multiselect("Select genes for this backbone", [g["gene_name"] for g in genes], key="design_selected_genes")
        option_cols = st.columns([2, 2, 3])
        with option_cols[0]:
            map_seq_type = st.selectbox(
                "Sequence type for this mapping",
                ["CDS", "Promoter", "Other"],
                index=["CDS", "Promoter", "Other"].index(default_seq_type),
                key="design_map_seq_type",
            )
        with option_cols[1]:
            map_remove_stop = st.checkbox(
                "Remove stop",
                value=bool(default_remove_stop),
                key="design_map_remove_stop",
                disabled=(map_seq_type != "CDS"),
                help="Only applies to CDS mappings.",
            )
        with option_cols[2]:
            st.caption("Add mapping with its own sequence type and stop-codon handling.")

        c1, c2, c3 = st.columns([1, 1, 2])
        with c1:
            if st.button("Add to batch", key="design_add_batch"):
                if not selected_genes:
                    st.warning("Please select at least one gene.")
                else:
                    effective_stop, stop_override = effective_remove_stop(selected_config, map_seq_type, map_remove_stop)
                    st.session_state["mapping_cart"].append(
                        {
                            **selected_source_ref,
                            "gene_names": selected_genes,
                            "seq_type": map_seq_type,
                            "remove_stop": effective_stop,
                            "remove_stop_override": stop_override,
                        }
                    )
                    st.success(
                        f"Added: {picked_label} -> {len(selected_genes)} genes "
                        f"[{map_seq_type}, remove_stop={_remove_stop_label(map_seq_type, effective_stop)}]"
                    )
        with c2:
            if st.button("Clear batch", key="design_clear_batch"):
                st.session_state["mapping_cart"] = []
                st.success("Cleared.")
        with c3:
            st.caption("Mapping mode keeps per-backbone parameter control.")

        cart_items = st.session_state["mapping_cart"]
    else:
        st.markdown("#### Cartesian build (Gene × Backbone)")
        bb_options = {r["backbone_label"]: r for r in source_refs}
        cart_backbone_labels = st.multiselect("Select backbones", list(bb_options.keys()), key="design_cart_backbones")
        cart_gene_names = st.multiselect("Select genes", [g["gene_name"] for g in genes], key="design_cart_genes")
        c1, c2 = st.columns(2)
        with c1:
            cart_seq_type = st.selectbox("Sequence type", ["CDS", "Promoter", "Other"], key="design_cart_seq_type")
        with c2:
            cart_remove_stop = st.checkbox(
                "Remove stop (CDS only)",
                value=bool(default_remove_stop),
                disabled=(cart_seq_type != "CDS"),
                key="design_cart_remove_stop",
            )
        cart_refs = [bb_options[x] for x in cart_backbone_labels]
        cart_items = _cartesian_items(cart_refs, cart_gene_names, cart_seq_type, cart_remove_stop)
        st.caption(f"Cartesian combinations: {len(cart_gene_names)} genes × {len(cart_refs)} backbones.")

    if cart_items:
        st.markdown("##### Current batch preview")
        preview_rows = _build_preview_rows(cart_items)
        st.dataframe(preview_rows, use_container_width=True)
    else:
        st.info("Batch is empty. Add at least one mapping/combination.")

    note = st.text_input("Run note (optional)", key="design_run_note")

    if st.button("Generate primers for current batch", key="design_run"):
        if not cart_items:
            st.warning("Batch is empty.")
        else:
            combos = _flatten_combos(cart_items)
            if not combos:
                st.warning("No valid combinations generated.")
            else:
                run_note = (note or "").strip()
                backbones_by_id = {int(b["id"]): b for b in backbones}
                try:
                    combo_configs = [_resolve_source_config(combo, backbones_by_id) for combo in combos]
                except (BackbonePresetError, ValueError) as exc:
                    st.error(f"Generate blocked: selected backbone configuration is invalid: {exc}")
                    conn.close()
                    return
                persist_run = all(c["backbone_source"] == "user_backbone" for c in combo_configs)
                params = dict(
                    batch_mode=batch_mode,
                    combo_count=len(combos),
                    min_len=int(min_len),
                    max_len=int(max_len),
                    tm_min=float(tm_min),
                    tm_max=float(tm_max),
                    target_tm=float(target_tm),
                    Na=float(Na),
                    Mg=float(Mg),
                    dNTPs=float(dNTPs),
                    primer_conc=float(primer_conc),
                    tm_calc_method="Biopython_Tm_NN",
                    enable_frame_check=bool(enable_frame_check),
                    construct_type=construct_type,
                    mappings=cart_items,
                    backbone_provenance=[config_provenance(c, combo["remove_stop"]) for c, combo in zip(combo_configs, combos)],
                    **policy_defaults,
                )
                run_id = None
                if persist_run:
                    first_backbone_id = int(combos[0]["backbone_id"])
                    cur = conn.cursor()
                    cur.execute(
                    """
                    INSERT INTO primer_runs(user_id, backbone_id, params_json, created_at, note, status)
                    VALUES(?,?,?,?,?,?)
                    """,
                        (
                            user["id"],
                            first_backbone_id,
                            json.dumps(params, ensure_ascii=False),
                            time.strftime("%Y-%m-%d %H:%M:%S"),
                            run_note or None,
                            "pending",
                        ),
                    )
                    run_id = int(cur.lastrowid)
                    conn.commit()
                else:
                    st.info("Shared preset results remain in this session; no shared preset is inserted into the user backbone database.")

                all_results = []
                failures = []
                all_genes = conn.execute(
                    "SELECT id, gene_name, sequence FROM genes WHERE user_id=? ORDER BY id",
                    (user["id"],),
                ).fetchall()
                gene_by_name = {g["gene_name"]: g for g in all_genes}
                total_steps = max(1, len(combos))
                done = 0
                prog = st.progress(0.0)

                for combo in combos:
                    bb = backbones_by_id.get(int(combo["backbone_id"])) if combo.get("backbone_id") is not None else None
                    try:
                        design_config = _resolve_source_config(combo, backbones_by_id)
                    except (BackbonePresetError, ValueError) as exc:
                        failures.append(
                            {
                                "backbone": combo.get("backbone_label", combo.get("preset_id", "unknown")),
                                "gene": combo["gene_name"],
                                "seq_type": combo["seq_type"],
                                "remove_stop": _remove_stop_label(combo["seq_type"], bool(combo["remove_stop"])),
                                "reason": str(exc),
                            }
                        )
                        done += 1
                        prog.progress(done / total_steps)
                        continue
                    g = gene_by_name.get(combo["gene_name"])
                    if bb is None and design_config["backbone_source"] == "user_backbone":
                        failures.append(
                            {
                                "backbone": combo.get("backbone_label", str(combo.get("backbone_id"))),
                                "gene": combo["gene_name"],
                                "seq_type": combo["seq_type"],
                                "remove_stop": _remove_stop_label(combo["seq_type"], bool(combo["remove_stop"])),
                                "reason": "Backbone no longer exists.",
                            }
                        )
                        done += 1
                        prog.progress(done / total_steps)
                        continue
                    if g is None:
                        failures.append(
                            {
                                "backbone": combo["backbone_label"],
                                "gene": combo["gene_name"],
                                "seq_type": combo["seq_type"],
                                "remove_stop": _remove_stop_label(combo["seq_type"], bool(combo["remove_stop"])),
                                "reason": "Gene no longer exists in database.",
                            }
                        )
                        done += 1
                        prog.progress(done / total_steps)
                        continue

                    try:
                        out = design_primers_for_gene(
                            gene_name=g["gene_name"],
                            seq=g["sequence"],
                            homo_f=design_config["homo_f"],
                            homo_r=design_config["homo_r"],
                            seq_type=combo["seq_type"],
                            remove_stop=bool(combo["remove_stop"]),
                            min_len=int(min_len),
                            max_len=int(max_len),
                            tm_min=float(tm_min),
                            tm_max=float(tm_max),
                            target_tm=float(target_tm),
                            gc_min=float(policy_defaults["gc_min"]),
                            gc_max=float(policy_defaults["gc_max"]),
                            Na=float(Na),
                            Mg=float(Mg),
                            dNTPs=float(dNTPs),
                            primer_conc=float(primer_conc),
                            length_mode=str(policy_defaults["length_mode"]),
                            cds_len_pref=tuple(policy_defaults["cds_len_pref"]),
                            non_cds_len_pref=tuple(policy_defaults["non_cds_len_pref"]),
                            tm_priority_weight=float(policy_defaults["tm_priority_weight"]),
                            non_cds_forbid_homopolymer=bool(policy_defaults["non_cds_forbid_homopolymer"]),
                            tm_pair_delta_max=float(policy_defaults["tm_pair_delta_max"]),
                            tm_delta_mode=str(policy_defaults["tm_delta_mode"]),
                            tm_profile=str(policy_defaults["tm_profile"]),
                            cds_recovery_enable=bool(policy_defaults["cds_recovery_enable"]),
                            enable_frame_check=(
                                bool(design_config["primer_engine_frame_check"])
                                if design_config["backbone_source"] == "shared_preset"
                                else bool(enable_frame_check)
                            ),
                            construct_type=design_config["construct_type"] if design_config["backbone_source"] == "shared_preset" else construct_type,
                            restriction_site_f=design_config["restriction_site_f"],
                            restriction_site_r=design_config["restriction_site_r"],
                            junction_extra_bases_f=int(design_config["junction_extra_bases_f"] or 0),
                            junction_extra_bases_r=int(design_config["junction_extra_bases_r"] or 0),
                        )
                        construct_name = build_construct_name(design_config["backbone_name"], g["gene_name"])
                        out["construct_name"] = construct_name
                        out["plasmid_name"] = construct_name
                        out["output_filename"] = build_map_filename(construct_name)
                        out["primer_name_f"] = build_primer_label(g["gene_name"], "F", backbone=design_config["backbone_name"])
                        out["primer_name_r"] = build_primer_label(g["gene_name"], "R", backbone=design_config["backbone_name"])
                        out["gene_id"] = g["id"]
                        out["backbone_id"] = design_config.get("backbone_id")
                        out["backbone_name"] = design_config["backbone_name"]
                        out["backbone_label"] = design_config["backbone_label"]
                        out["gene_name"] = g["gene_name"]
                        out["seq_type"] = combo["seq_type"]
                        out["remove_stop"] = bool(combo["remove_stop"])
                        out["construct_type"] = design_config["construct_type"] if design_config["backbone_source"] == "shared_preset" else construct_type
                        out["backbone_source"] = design_config["backbone_source"]
                        out["preset_id"] = design_config.get("preset_id")
                        out["preset_version"] = design_config.get("preset_version")
                        out["vector_sha256"] = design_config.get("vector_sha256")
                        out["fusion_order"] = design_config.get("fusion_order")
                        out["fusion_architecture"] = design_config.get("fusion_architecture")
                        out["validated_homology_arm_length"] = design_config.get("homology_arm_length")
                        out["fusion_frame_validation"] = design_config.get("fusion_frame_validation")
                        out["homo_f"] = design_config["homo_f"]
                        out["homo_r"] = design_config["homo_r"]
                        out["vector_sequence"] = design_config["vector_sequence"]
                        out["insertion_mode"] = design_config["insertion_mode"]
                        out["restriction_site_f"] = design_config["restriction_site_f"]
                        out["restriction_site_r"] = design_config["restriction_site_r"]
                        out["junction_extra_bases_f"] = int(design_config["junction_extra_bases_f"] or 0)
                        out["junction_extra_bases_r"] = int(design_config["junction_extra_bases_r"] or 0)
                        out["vector_map_note"] = design_config["vector_map_note"]
                        out["vector_features_json"] = design_config["vector_features_json"]
                        out["vector_topology"] = design_config["vector_topology"]
                        out["vector_source_format"] = design_config["vector_source_format"]
                        out["vector_record_name"] = design_config["vector_record_name"]
                        out["backbone_provenance"] = config_provenance(design_config, bool(combo["remove_stop"]))
                        all_results.append(out)

                        if run_id is not None:
                            conn.execute(
                            """
                            INSERT INTO primers(
                                run_id, gene_id, backbone_id,
                                primer_name_f, primer_name_r,
                                gene_bind_f, gene_bind_r, final_f, final_r,
                                bind_start_f, bind_end_f, bind_start_r, bind_end_r,
                                amplicon_seq, insert_seq, seq_type, remove_stop, construct_type,
                                tm_f, tm_r, gc_f, gc_r, amplicon_len, warnings, selected_reason,
                                tm_delta, rescue_tier, tm_profile, tm_calc_params
                            )
                            VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                            """,
                            (
                                run_id,
                                g["id"],
                                int(bb["id"]),
                                out["primer_name_f"],
                                out["primer_name_r"],
                                out["gene_bind_f"],
                                out["gene_bind_r"],
                                out["final_f"],
                                out["final_r"],
                                out.get("bind_start_f"),
                                out.get("bind_end_f"),
                                out.get("bind_start_r"),
                                out.get("bind_end_r"),
                                out.get("amplicon_seq", ""),
                                out.get("insert_seq", ""),
                                out.get("seq_type", ""),
                                1 if bool(out.get("remove_stop")) else 0,
                                out.get("construct_type", ""),
                                out.get("tm_f"),
                                out.get("tm_r"),
                                out.get("gc_f"),
                                out.get("gc_r"),
                                out.get("amplicon_len"),
                                out.get("warnings", ""),
                                out.get("selected_reason", ""),
                                out.get("tm_delta"),
                                out.get("rescue_tier", ""),
                                out.get("tm_profile", ""),
                                json.dumps(out.get("tm_calc_params", {}), ensure_ascii=False),
                            ),
                            )
                            conn.commit()
                    except Exception as exc:
                        failures.append(
                            {
                                "backbone": design_config.get("backbone_label", "unknown"),
                                "gene": g["gene_name"],
                                "seq_type": combo["seq_type"],
                                "remove_stop": _remove_stop_label(combo["seq_type"], bool(combo["remove_stop"])),
                                "reason": str(exc),
                            }
                        )
                    done += 1
                    prog.progress(done / total_steps)

                run_status = "success" if not failures else "fail"
                if run_id is not None:
                    conn.execute("UPDATE primer_runs SET status=? WHERE id=? AND user_id=?", (run_status, run_id, user["id"]))
                    conn.commit()

                st.session_state["last_results"] = all_results
                st.session_state["last_failures"] = failures
                st.session_state["last_run_id"] = run_id
                run_label = str(run_id) if run_id is not None else "session-only shared preset run"
                st.success(f"Generated primers. Run {run_label} status: {run_status}.")
                if failures:
                    st.warning(f"{len(failures)} combinations failed and were skipped.")

    if "last_results" in st.session_state and st.session_state["last_results"]:
        st.divider()
        st.markdown("### Latest results (this session)")
        show_warning_only = st.checkbox("Show warning-only rows", key="design_warning_only")
        all_result_rows = st.session_state["last_results"]
        result_rows = all_result_rows
        if show_warning_only:
            result_rows = [r for r in result_rows if (r.get("warnings") or "").strip()]
        st.dataframe(
            [
                {
                    k: v
                    for k, v in r.items()
                    if k
                    in [
                        "backbone_label",
                        "gene_name",
                        "primer_name_f",
                        "primer_name_r",
                        "seq_type",
                        "construct_type",
                        "remove_stop",
                        "amplicon_len",
                        "tm_f",
                        "tm_r",
                        "tm_delta",
                        "gc_f",
                        "gc_r",
                        "tm_profile",
                        "rescue_tier",
                        "selected_reason",
                        "warnings",
                    ]
                }
                for r in result_rows
            ],
            use_container_width=True,
        )
        try:
            validate_primer_rows(st.session_state["last_results"])
            template_rows = get_user_templates(conn, int(user["id"]))
            options = [("built_in", "Built-in template")]
            default_idx = 0
            default_tid = get_default_template_id(conn, int(user["id"]))
            for r in template_rows:
                tid = str(int(r["id"]))
                label = f"{r['template_name']} (#{r['id']})"
                if int(r["is_default"]) == 1:
                    label += " [default]"
                options.append((tid, label))
                if default_tid is not None and int(r["id"]) == int(default_tid):
                    default_idx = len(options) - 1
            selected_template_key = st.selectbox(
                "Excel template",
                [x[0] for x in options],
                index=default_idx,
                format_func=lambda k: dict(options).get(k, k),
                key="design_export_template_select",
            )
            if selected_template_key == "built_in":
                resolved = resolve_order_template(template_path=GLOBAL_DEFAULT_TEMPLATE_PATH)
            else:
                resolved = resolve_order_template(
                    template_id=int(selected_template_key),
                    user_id=int(user["id"]),
                    conn=conn,
                )
            xbytes = export_from_template_bytes(
                st.session_state["last_results"],
                template_bytes=resolved.template_bytes,
                profile=profile,
            )
            st.download_button(
                "Download Sangon Excel (Template)",
                data=xbytes,
                file_name="Sangon_Primers_Template.xlsx",
                key="design_download_excel",
            )
            st.caption(f"Export template: {resolved.template_name}")
        except Exception as exc:
            st.error(f"Excel export blocked by validation: {exc}")

        frag_txt = build_fragment_length_txt(all_result_rows)
        st.download_button(
            "Download Fragment Length TXT",
            data=frag_txt,
            file_name="Fragment_Lengths.txt",
            key="design_download_frag_txt",
        )
        frag_simple_txt = build_fragment_length_simple_txt(all_result_rows)
        st.download_button(
            "Download Fragment Length TXT (Simple)",
            data=frag_simple_txt,
            file_name="Fragment_Lengths_Simple.txt",
            key="design_download_frag_simple_txt",
        )
        non_cds_fasta = build_non_cds_amplicon_fasta(all_result_rows)
        if non_cds_fasta:
            st.download_button(
                "Download Non-CDS Amplicons FASTA",
                data=non_cds_fasta,
                file_name="NonCDS_Amplicons.fasta",
                key="design_download_non_cds_fasta",
            )
        else:
            st.info("No non-CDS amplicon sequences in this batch.")

        sequencing_download(all_result_rows, "design_sequencing")

        st.markdown("#### Construct map export")
        first_seq_type = str(all_result_rows[0].get("seq_type") or "CDS").lower()
        default_feature = "CDS" if first_seq_type == "cds" else ("promoter" if first_seq_type == "promoter" else "misc_feature")
        m1, m2, m3 = st.columns([1, 2, 1])
        with m1:
            insert_feature_type = st.selectbox(
                "Insert feature type",
                ["CDS", "promoter", "misc_feature"],
                index=["CDS", "promoter", "misc_feature"].index(default_feature),
                key="design_map_insert_feature_type",
            )
        with m2:
            default_feature_name = all_result_rows[0].get("gene_name", "") if len(all_result_rows) == 1 else ""
            insert_feature_name = st.text_input(
                "Insert feature name",
                value=default_feature_name,
                placeholder="Blank = use each gene name",
                key="design_map_insert_feature_name",
            )
        with m3:
            export_genbank = st.checkbox("Export GenBank map", value=True, key="design_map_export_genbank")
        with st.expander("Advanced map options", expanded=False):
            insert_feature_note = st.text_input("Insert feature note/color", value="", key="design_map_insert_feature_note")
            homology_feature_type = st.selectbox(
                "Homology arm feature type",
                ["misc_feature", "primer_bind"],
                key="design_map_homology_feature_type",
            )
            st.caption("Insert direction is inferred from the homology arms; if inference fails, the export will show a Chinese action hint.")

        map_rows = _with_map_options(all_result_rows, insert_feature_type, insert_feature_name, insert_feature_note, homology_feature_type)
        if export_genbank:
            try:
                gb_name, gb_data, gb_mime = build_genbank_bundle(map_rows)
                st.download_button("Download GenBank map", data=gb_data, file_name=gb_name, mime=gb_mime, key="design_download_genbank")
            except Exception as exc:
                st.warning(f"GenBank 图谱导出失败：{exc}")

        snapgene_exe, _ = find_snapgene_cli()
        if snapgene_exe is None:
            st.info("当前环境未检测到 SnapGene Desktop CLI，已提供 GenBank `.gbk`；本地安装 SnapGene 后可导出 `.dna`。")
        else:
            try:
                dna_name, dna_data, dna_mime, used_exe = build_snapgene_dna_bundle(map_rows)
                st.download_button(
                    "Download SnapGene DNA map",
                    data=dna_data,
                    file_name=dna_name,
                    mime=dna_mime,
                    key="design_download_snapgene_dna",
                )
                st.caption(f"SnapGene CLI: {used_exe}")
            except Exception as exc:
                st.warning(f"SnapGene `.dna` 转换失败：{exc}")

    failures = st.session_state.get("last_failures", [])
    if failures:
        st.markdown("### Failed combinations (this session)")
        st.dataframe(failures, use_container_width=True)
        failed_gene_list = "\n".join([f"{f['backbone']}::{f['gene']}" for f in failures])
        st.text_area("Failed items (copy for re-run)", value=failed_gene_list, height=120, key="design_failed_copy")

    conn.close()
