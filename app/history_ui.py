from sequencing_ui import sequencing_download
import json
import math

import streamlit as st

from db import get_conn
from excel_export import export_from_template_bytes, validate_primer_rows
from fragment_export import build_fragment_length_simple_txt, build_fragment_length_txt, build_non_cds_amplicon_fasta
from genbank_export import build_genbank_bundle
from snapgene_export import build_snapgene_dna_bundle, find_snapgene_cli
from templates_ui import get_default_template_id, get_effective_template_bytes, get_user_templates
from template_resolver import GLOBAL_DEFAULT_TEMPLATE_PATH, resolve_order_template
from naming import build_map_filename, construct_from_primer_name


def _fmt_tm(v, label):
    if v is None:
        return f"{label} NA"
    try:
        if v != v:
            return f"{label} NA"
    except Exception:
        return f"{label} NA"
    try:
        vf = float(v)
        if math.isnan(vf):
            return f"{label} NA"
        return f"{label} {vf:.2f}"
    except Exception:
        return f"{label} NA"


def _fmt_delta(v):
    if v is None:
        return "dTm NA"
    try:
        if v != v:
            return "dTm NA"
    except Exception:
        return "dTm NA"
    try:
        vf = float(v)
        if math.isnan(vf):
            return "dTm NA"
        return f"dTm {vf:.2f}"
    except Exception:
        return "dTm NA"


def _bb_label(r) -> str:
    enz = (r["enzyme_sites"] or "").strip()
    return f"{r['backbone_name']} [{enz}]" if enz else f"{r['backbone_name']} [no-enzyme]"


def _with_map_options(primer_rows, feature_type: str, feature_name: str, feature_note: str, homology_type: str):
    out = []
    for row in primer_rows:
        item = dict(row)
        item["insert_feature_type"] = feature_type
        item["insert_feature_name"] = (feature_name or "").strip() or item.get("gene_name") or "insert"
        item["insert_feature_note"] = (feature_note or "").strip() or "Inserted sequence"
        item["homology_feature_type"] = homology_type
        out.append(item)
    return out


def history_page(user, config):
    st.subheader("History")
    conn = get_conn()

    runs = conn.execute(
        """
        SELECT r.id, r.created_at, r.note, r.params_json, r.status,
               b.name AS backbone_name, b.enzyme_sites
        FROM primer_runs r
        LEFT JOIN backbones b ON b.id = r.backbone_id
        WHERE r.user_id=?
        ORDER BY r.id DESC
        LIMIT 200
        """,
        (user["id"],),
    ).fetchall()

    if not runs:
        st.info("No runs yet.")
        conn.close()
        return

    run_ids = [r["id"] for r in runs]
    selected_run = st.selectbox("Select a run", run_ids, format_func=lambda rid: f"Run {rid}", key="history_run_pick")
    rrow = next(r for r in runs if r["id"] == selected_run)

    status = str(rrow["status"] or "pending")
    st.write(f"Status: **{status}**")
    if rrow["backbone_name"]:
        st.write(f"Backbone: **{_bb_label(rrow)}**")
    else:
        st.write("Backbone: **mixed/legacy**")
    st.write(f"Created: {rrow['created_at']}")
    if rrow["note"]:
        st.write(f"Note: {rrow['note']}")
    try:
        params = json.loads(rrow["params_json"] or "{}")
    except Exception:
        params = {}
    if params.get("batch_mode"):
        st.caption(f"Batch mode: {params.get('batch_mode')}")
    if params.get("combo_count"):
        st.caption(f"Combinations: {params.get('combo_count')}")
    if params.get("tm_profile"):
        st.caption(f"Tm profile: {params.get('tm_profile')}")
    if params.get("construct_type"):
        st.caption(f"Construct type: {params.get('construct_type')}")

    dcol1, dcol2, dcol3 = st.columns([2, 3, 3])
    with dcol1:
        if st.button("Delete this run", key="history_delete_current"):
            conn.execute("DELETE FROM primer_runs WHERE user_id=? AND id=?", (user["id"], selected_run))
            conn.commit()
            st.success(f"Deleted run {selected_run}.")
            st.rerun()
    with dcol2:
        batch_del = st.multiselect("Batch delete runs", run_ids, key="history_delete_batch_pick")
    with dcol3:
        confirm_batch = st.checkbox("Confirm batch deletion", key="history_delete_batch_confirm")
        if st.button("Delete selected runs", key="history_delete_batch_btn"):
            if not batch_del:
                st.warning("No runs selected.")
            elif not confirm_batch:
                st.warning("Please confirm batch deletion first.")
            else:
                placeholders = ",".join("?" for _ in batch_del)
                conn.execute(f"DELETE FROM primer_runs WHERE user_id=? AND id IN ({placeholders})", [user["id"], *batch_del])
                conn.commit()
                st.success(f"Deleted {len(batch_del)} run(s).")
                st.rerun()

    prof = conn.execute("SELECT * FROM user_profiles WHERE user_id=?", (user["id"],)).fetchone()
    profile = dict(prof) if prof else None

    primers = conn.execute(
        """
        SELECT
            p.*,
            g.gene_name,
            b.name AS backbone_name,
            b.enzyme_sites,
            b.homo_f,
            b.homo_r,
            b.vector_sequence,
            b.insertion_mode,
            b.restriction_site_f,
            b.restriction_site_r,
            b.junction_extra_bases_f,
            b.junction_extra_bases_r,
            b.vector_map_note,
            b.vector_features_json,
            b.vector_topology,
            b.vector_source_format,
            b.vector_record_name
        FROM primers p
        JOIN genes g ON g.id = p.gene_id
        JOIN primer_runs r ON r.id = p.run_id
        LEFT JOIN backbones b ON b.id = p.backbone_id
        WHERE p.run_id=? AND r.user_id=? AND g.user_id=?
        ORDER BY p.id ASC
        """,
        (selected_run, user["id"], user["id"]),
    ).fetchall()

    st.markdown("#### Primers")
    for p in primers[:200]:
        cols = st.columns([3, 3, 2, 2, 2, 3])
        cols[0].write(f"**{p['gene_name']}**")
        bb_name = p["backbone_name"] if p["backbone_name"] else "N/A"
        cols[1].write(bb_name)
        cols[2].write(_fmt_tm(p["tm_f"], "TmF"))
        cols[3].write(_fmt_tm(p["tm_r"], "TmR"))
        cols[4].write(_fmt_delta(p["tm_delta"]))
        cols[5].write(f"Amp {p['amplicon_len'] if p['amplicon_len'] is not None else 'NA'} bp")
        st.code(f"{p['primer_name_f']}: {p['final_f']}\n{p['primer_name_r']}: {p['final_r']}", language="text")
        if p["warnings"]:
            st.caption(f"Warnings: {p['warnings']}")

    if len(primers) > 200:
        st.info(f"Showing first 200 primer pairs out of {len(primers)} pairs.")

    primer_rows = []
    for p in primers:
        construct_name = construct_from_primer_name(p["primer_name_f"])
        primer_rows.append(
            {
                "primer_name_f": p["primer_name_f"],
                "primer_name_r": p["primer_name_r"],
                "final_f": p["final_f"],
                "final_r": p["final_r"],
                "gene_bind_f": p["gene_bind_f"],
                "gene_bind_r": p["gene_bind_r"],
                "bind_start_f": p["bind_start_f"],
                "bind_end_f": p["bind_end_f"],
                "bind_start_r": p["bind_start_r"],
                "bind_end_r": p["bind_end_r"],
                "amplicon_seq": p["amplicon_seq"] or "",
                "insert_seq": p["insert_seq"] or "",
                "backbone_name": p["backbone_name"] or "",
                "gene_name": p["gene_name"],
                "construct_name": construct_name,
                "plasmid_name": construct_name,
                "output_filename": build_map_filename(construct_name),
                "seq_type": p["seq_type"] or params.get("seq_type", ""),
                "remove_stop": bool(p["remove_stop"]) if p["remove_stop"] is not None else bool(params.get("remove_stop", False)),
                "construct_type": p["construct_type"] or params.get("construct_type", ""),
                "amplicon_len": p["amplicon_len"],
                "tm_f": p["tm_f"],
                "tm_r": p["tm_r"],
                "tm_delta": p["tm_delta"],
                "tm_profile": p["tm_profile"] or params.get("tm_profile", ""),
                "rescue_tier": p["rescue_tier"] or "",
                "warnings": p["warnings"] or "",
                "selected_reason": p["selected_reason"] or "",
                "homo_f": p["homo_f"] or "",
                "homo_r": p["homo_r"] or "",
                "vector_sequence": p["vector_sequence"] or "",
                "insertion_mode": p["insertion_mode"] or "linear_join",
                "restriction_site_f": p["restriction_site_f"] or "",
                "restriction_site_r": p["restriction_site_r"] or "",
                "junction_extra_bases_f": int(p["junction_extra_bases_f"] or 0),
                "junction_extra_bases_r": int(p["junction_extra_bases_r"] or 0),
                "vector_map_note": p["vector_map_note"] or "",
                "vector_features_json": p["vector_features_json"] or "[]",
                "vector_topology": p["vector_topology"] or "",
                "vector_source_format": p["vector_source_format"] or "",
                "vector_record_name": p["vector_record_name"] or "",
            }
        )

    try:
        validate_primer_rows(primer_rows)
        st.success("Re-export consistency check passed.")
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
            "Excel template for re-export",
            [x[0] for x in options],
            index=default_idx,
            format_func=lambda k: dict(options).get(k, k),
            key=f"history_export_template_{selected_run}",
        )
        if selected_template_key == "built_in":
            resolved = resolve_order_template(template_path=GLOBAL_DEFAULT_TEMPLATE_PATH)
        else:
            resolved = resolve_order_template(
                template_id=int(selected_template_key),
                user_id=int(user["id"]),
                conn=conn,
            )
        xbytes = export_from_template_bytes(primer_rows, template_bytes=resolved.template_bytes, profile=profile)
        st.download_button(
            "Download Sangon Excel (Template)",
            data=xbytes,
            file_name=f"Sangon_Run_{selected_run}.xlsx",
            key=f"history_download_{selected_run}",
        )
        st.caption(f"Export template: {resolved.template_name}")
        frag_txt = build_fragment_length_txt(primer_rows)
        st.download_button(
            "Download Fragment Length TXT",
            data=frag_txt,
            file_name=f"Fragment_Lengths_Run_{selected_run}.txt",
            key=f"history_download_frag_txt_{selected_run}",
        )
        frag_simple_txt = build_fragment_length_simple_txt(primer_rows)
        st.download_button(
            "Download Fragment Length TXT (Simple)",
            data=frag_simple_txt,
            file_name=f"Fragment_Lengths_Simple_Run_{selected_run}.txt",
            key=f"history_download_frag_simple_txt_{selected_run}",
        )
        non_cds_fasta = build_non_cds_amplicon_fasta(primer_rows)
        if non_cds_fasta:
            st.download_button(
                "Download Non-CDS Amplicons FASTA",
                data=non_cds_fasta,
                file_name=f"NonCDS_Amplicons_Run_{selected_run}.fasta",
                key=f"history_download_non_cds_fasta_{selected_run}",
            )
        else:
            st.info("This run has no non-CDS amplicon sequences.")

        sequencing_download(primer_rows, f"history_sequencing_{selected_run}")

        st.markdown("#### Construct map export")
        first_seq_type = str(primer_rows[0].get("seq_type") or "CDS").lower() if primer_rows else "cds"
        default_feature = "CDS" if first_seq_type == "cds" else ("promoter" if first_seq_type == "promoter" else "misc_feature")
        h1, h2, h3 = st.columns([1, 2, 1])
        with h1:
            insert_feature_type = st.selectbox(
                "Insert feature type",
                ["CDS", "promoter", "misc_feature"],
                index=["CDS", "promoter", "misc_feature"].index(default_feature),
                key=f"history_map_insert_feature_type_{selected_run}",
            )
        with h2:
            default_feature_name = primer_rows[0].get("gene_name", "") if len(primer_rows) == 1 else ""
            insert_feature_name = st.text_input(
                "Insert feature name",
                value=default_feature_name,
                placeholder="Blank = use each gene name",
                key=f"history_map_insert_feature_name_{selected_run}",
            )
        with h3:
            export_genbank = st.checkbox("Export GenBank map", value=True, key=f"history_map_export_genbank_{selected_run}")
        with st.expander("Advanced map options", expanded=False):
            insert_feature_note = st.text_input("Insert feature note/color", value="", key=f"history_map_insert_feature_note_{selected_run}")
            homology_feature_type = st.selectbox(
                "Homology arm feature type",
                ["misc_feature", "primer_bind"],
                key=f"history_map_homology_feature_type_{selected_run}",
            )

        map_rows = _with_map_options(primer_rows, insert_feature_type, insert_feature_name, insert_feature_note, homology_feature_type)
        if export_genbank:
            gb_name, gb_data, gb_mime = build_genbank_bundle(map_rows)
            st.download_button(
                "Download GenBank map",
                data=gb_data,
                file_name=f"Run_{selected_run}_{gb_name}",
                mime=gb_mime,
                key=f"history_download_genbank_{selected_run}",
            )

        snapgene_exe, _ = find_snapgene_cli()
        if snapgene_exe is None:
            st.info("当前环境未检测到 SnapGene Desktop CLI，已提供 GenBank `.gbk`；本地安装 SnapGene 后可导出 `.dna`。")
        else:
            try:
                dna_name, dna_data, dna_mime, used_exe = build_snapgene_dna_bundle(map_rows)
                st.download_button(
                    "Download SnapGene DNA map",
                    data=dna_data,
                    file_name=f"Run_{selected_run}_{dna_name}",
                    mime=dna_mime,
                    key=f"history_download_snapgene_dna_{selected_run}",
                )
                st.caption(f"SnapGene CLI: {used_exe}")
            except Exception as exc:
                st.warning(f"SnapGene `.dna` 转换失败：{exc}")
    except Exception as exc:
        st.error(f"Re-export blocked by validation: {exc}")

    conn.close()
