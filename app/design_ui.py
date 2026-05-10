import time
import json
import re
import streamlit as st
from db import get_conn
from primer_design import design_primers_for_gene
from excel_export import TEMPLATE_PATH, export_from_template_bytes, validate_primer_rows
from fragment_export import build_fragment_length_txt, build_fragment_length_simple_txt, build_non_cds_amplicon_fasta
from templates_ui import get_user_templates, get_default_template_id, get_effective_template_bytes


def _sanitize_name(s: str) -> str:
    s = (s or "").strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_.-]", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def _primer_name(backbone_name: str, gene_name: str, direction: str) -> str:
    b = _sanitize_name(backbone_name)
    g = _sanitize_name(gene_name)
    d = "F" if direction.upper().startswith("F") else "R"
    return f"{b}_{g}-{d}"


def _remove_stop_label(seq_type: str, remove_stop: bool) -> str:
    if (seq_type or "").upper() != "CDS":
        return "N/A"
    return "Yes" if remove_stop else "No"


def _bb_label(backbone_row) -> str:
    enz = (backbone_row["enzyme_sites"] or "").strip()
    return f"{backbone_row['name']} [{enz}]" if enz else f"{backbone_row['name']} [no-enzyme]"


def _normalize_mapping_cart(cart, backbones):
    by_id = {int(b["id"]): b for b in backbones}
    by_name_first = {}
    for b in backbones:
        by_name_first.setdefault(b["name"], b)

    fixed = []
    for item in cart:
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


def design_page(user, config):
    st.subheader("Design (Mapping-based batch)")
    conn = get_conn()

    prof = conn.execute("SELECT * FROM user_profiles WHERE user_id=?", (user["id"],)).fetchone()
    profile = dict(prof) if prof else None

    backbones = conn.execute(
        "SELECT * FROM backbones WHERE user_id=? ORDER BY name, enzyme_sites, id",
        (user["id"],),
    ).fetchall()
    if not backbones:
        st.warning("Please add at least one backbone first.")
        conn.close()
        return

    defaults = json.loads((backbones[0]["defaults_json"] or "{}"))

    st.markdown("#### Parameters")
    col1, col2, col3 = st.columns(3)
    with col1:
        default_seq_type = st.selectbox(
            "Default sequence type for new mapping",
            ["CDS", "Promoter", "Other"],
            key="design_default_seq_type",
        )
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
    with col3:
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
        "gc_min": 40.0,
        "gc_max": 60.0,
        "length_mode": "adaptive_soft",
        "cds_len_pref": [18, 21],
        "non_cds_len_pref": [20, 26],
        "tm_priority_weight": 2.5,
        "non_cds_forbid_homopolymer": True,
        "tm_pair_delta_max": float(tm_pair_delta_max),
        "tm_delta_mode": "soft",
        "tm_profile": tm_profile,
        "cds_recovery_enable": True,
    }
    st.caption("Design policy: Tm-priority + GC 40-60 hard filter + adaptive soft length preference + soft |TmF-TmR| + CDS auto recovery.")

    st.markdown("#### Build mapping (Backbone -> Genes)")
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

    if "mapping_cart" not in st.session_state:
        st.session_state["mapping_cart"] = []
    st.session_state["mapping_cart"] = _normalize_mapping_cart(st.session_state["mapping_cart"], backbones)

    backbone_options = { _bb_label(b): int(b["id"]) for b in backbones }
    picked_label = st.selectbox("Select backbone (one at a time)", list(backbone_options.keys()), key="design_backbone_pick")
    picked_id = backbone_options[picked_label]
    picked_bb = next(b for b in backbones if int(b["id"]) == picked_id)

    selected_genes = st.multiselect(
        "Select genes for this backbone",
        [g["gene_name"] for g in genes],
        key="design_selected_genes",
    )

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
                effective_remove_stop = bool(map_remove_stop) if map_seq_type == "CDS" else False
                st.session_state["mapping_cart"].append(
                    {
                        "backbone_id": int(picked_bb["id"]),
                        "backbone_name": picked_bb["name"],
                        "backbone_label": _bb_label(picked_bb),
                        "gene_names": selected_genes,
                        "seq_type": map_seq_type,
                        "remove_stop": effective_remove_stop,
                    }
                )
                st.success(
                    f"Added: {_bb_label(picked_bb)} -> {len(selected_genes)} genes "
                    f"[{map_seq_type}, remove_stop={_remove_stop_label(map_seq_type, effective_remove_stop)}]"
                )
    with c2:
        if st.button("Clear batch", key="design_clear_batch"):
            st.session_state["mapping_cart"] = []
            st.success("Cleared.")
    with c3:
        st.caption("You control pairing directly. No forced Cartesian product.")

    if st.session_state["mapping_cart"]:
        st.markdown("##### Current batch mapping")
        for idx, item in enumerate(st.session_state["mapping_cart"], start=1):
            cols = st.columns([3, 5, 3, 1])
            cols[0].write(f"{idx}. **{item['backbone_label']}**")
            cols[1].write(", ".join(item["gene_names"]))
            cols[2].write(
                f"{item['seq_type']} | remove_stop={_remove_stop_label(item['seq_type'], bool(item.get('remove_stop')))}"
            )
            if cols[3].button("X", key=f"rm_{idx}"):
                st.session_state["mapping_cart"].pop(idx - 1)
                st.rerun()
    else:
        st.info("Batch is empty. Add at least one mapping above.")

    preview_rows = _build_preview_rows(st.session_state["mapping_cart"])
    if preview_rows:
        st.markdown("##### Batch preview")
        st.dataframe(preview_rows, use_container_width=True)

    note = st.text_input("Run note (optional)", key="design_run_note")

    if st.button("Generate primers for batch mapping", key="design_run"):
        if not st.session_state["mapping_cart"]:
            st.warning("Batch is empty.")
        else:
            group_tag = time.strftime("%Y%m%d_%H%M%S")
            base_note = (note or "").strip()
            run_ids = []
            all_results = []
            failures = []

            all_genes = conn.execute(
                "SELECT id, gene_name, sequence FROM genes WHERE user_id=? ORDER BY id",
                (user["id"],),
            ).fetchall()
            gene_by_name = {g["gene_name"]: g for g in all_genes}

            total_steps = sum(len(it["gene_names"]) for it in st.session_state["mapping_cart"])
            total_steps = max(1, total_steps)
            done = 0
            prog = st.progress(0.0)

            backbones_by_id = {int(b["id"]): b for b in backbones}

            for item in st.session_state["mapping_cart"]:
                bb = backbones_by_id.get(int(item["backbone_id"]))
                if bb is None:
                    failures.append(
                        {
                            "backbone": item.get("backbone_label", str(item.get("backbone_id"))),
                            "gene": ",".join(item.get("gene_names", [])),
                            "seq_type": item.get("seq_type", ""),
                            "remove_stop": _remove_stop_label(item.get("seq_type", ""), bool(item.get("remove_stop"))),
                            "reason": "Backbone no longer exists.",
                        }
                    )
                    done += len(item.get("gene_names", []))
                    prog.progress(done / total_steps)
                    continue

                gene_list = [gene_by_name[n] for n in item["gene_names"] if n in gene_by_name]
                missing_names = [n for n in item["gene_names"] if n not in gene_by_name]
                for mn in missing_names:
                    failures.append(
                        {
                            "backbone": item["backbone_label"],
                            "gene": mn,
                            "seq_type": item.get("seq_type", ""),
                            "remove_stop": _remove_stop_label(item.get("seq_type", ""), bool(item.get("remove_stop"))),
                            "reason": "Gene no longer exists in database.",
                        }
                    )

                item_seq_type = str(item.get("seq_type") or "CDS")
                item_remove_stop = bool(item.get("remove_stop")) if item_seq_type.upper() == "CDS" else False

                params = dict(
                    seq_type=item_seq_type,
                    remove_stop=item_remove_stop,
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
                    **policy_defaults,
                )

                cur = conn.cursor()
                run_note = f"[group:{group_tag}] " + (base_note if base_note else "")
                cur.execute(
                    "INSERT INTO primer_runs(user_id, backbone_id, params_json, created_at, note) VALUES(?,?,?,?,?)",
                    (
                        user["id"],
                        bb["id"],
                        json.dumps(params, ensure_ascii=False),
                        time.strftime("%Y-%m-%d %H:%M:%S"),
                        run_note or None,
                    ),
                )
                run_id = cur.lastrowid
                conn.commit()
                run_ids.append(run_id)

                for g in gene_list:
                    try:
                        out = design_primers_for_gene(
                            gene_name=g["gene_name"],
                            seq=g["sequence"],
                            homo_f=bb["homo_f"],
                            homo_r=bb["homo_r"],
                            seq_type=item_seq_type,
                            remove_stop=item_remove_stop,
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
                        )
                        out["primer_name_f"] = _primer_name(bb["name"], g["gene_name"], "F")
                        out["primer_name_r"] = _primer_name(bb["name"], g["gene_name"], "R")
                        out["gene_id"] = g["id"]
                        out["backbone_name"] = bb["name"]
                        out["backbone_label"] = _bb_label(bb)
                        out["gene_name"] = g["gene_name"]
                        out["seq_type"] = item_seq_type
                        out["remove_stop"] = item_remove_stop
                        all_results.append(out)

                        conn.execute(
                            """
                            INSERT INTO primers(run_id, gene_id, primer_name_f, primer_name_r,
                                                gene_bind_f, gene_bind_r, final_f, final_r,
                                                bind_start_f, bind_end_f, bind_start_r, bind_end_r, amplicon_seq,
                                                tm_f, tm_r, gc_f, gc_r, amplicon_len, warnings, selected_reason,
                                                tm_delta, rescue_tier, tm_profile, tm_calc_params)
                            VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                            """,
                            (
                                run_id,
                                g["id"],
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
                    except Exception as e:
                        failures.append(
                            {
                                "backbone": _bb_label(bb),
                                "gene": g["gene_name"],
                                "seq_type": item_seq_type,
                                "remove_stop": _remove_stop_label(item_seq_type, item_remove_stop),
                                "reason": str(e),
                            }
                        )
                    done += 1
                    prog.progress(done / total_steps)

            st.session_state["last_results"] = all_results
            st.session_state["last_failures"] = failures
            st.success(f"Generated primers for mapping. Saved runs: {run_ids}")
            if failures:
                st.warning(f"{len(failures)} genes failed and were skipped.")

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
                    if k in [
                        "backbone_label",
                        "gene_name",
                        "primer_name_f",
                        "primer_name_r",
                        "seq_type",
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
                template_bytes = TEMPLATE_PATH.read_bytes()
                template_label = "Built-in template"
            else:
                template_bytes, template_label = get_effective_template_bytes(
                    conn, int(user["id"]), int(selected_template_key)
                )
            xbytes = export_from_template_bytes(
                st.session_state["last_results"],
                template_bytes=template_bytes,
                profile=profile,
            )
            st.download_button(
                "Download Sangon Excel (Template)",
                data=xbytes,
                file_name="Sangon_Primers_Template.xlsx",
                key="design_download_excel",
            )
            st.caption(f"Export template: {template_label}")
        except Exception as e:
            st.error(f"Excel export blocked by validation: {e}")
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

    failures = st.session_state.get("last_failures", [])
    if failures:
        st.markdown("### Failed genes (this session)")
        st.dataframe(failures, use_container_width=True)
        failed_gene_list = "\n".join([f["gene"] for f in failures])
        st.text_area("Failed gene names (copy for re-run)", value=failed_gene_list, height=120, key="design_failed_copy")

    conn.close()

