import json
import streamlit as st
from db import get_conn
from excel_export import TEMPLATE_PATH, export_from_template_bytes, validate_primer_rows
from fragment_export import build_fragment_length_txt, build_fragment_length_simple_txt, build_non_cds_amplicon_fasta
from templates_ui import get_user_templates, get_default_template_id, get_effective_template_bytes


def _fmt_tm(v, label):
    if v is None:
        return f"{label} NA"
    try:
        if v != v:
            return f"{label} NA"
    except Exception:
        return f"{label} NA"
    try:
        return f"{label} {float(v):.2f}"
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
        return f"dTm {float(v):.2f}"
    except Exception:
        return "dTm NA"


def _bb_label(r) -> str:
    enz = (r["enzyme_sites"] or "").strip()
    return f"{r['backbone_name']} [{enz}]" if enz else f"{r['backbone_name']} [no-enzyme]"


def history_page(user, config):
    st.subheader("History")
    conn = get_conn()

    runs = conn.execute(
        """
        SELECT r.id, r.created_at, r.note, r.params_json, b.name AS backbone_name, b.enzyme_sites
        FROM primer_runs r
        JOIN backbones b ON b.id = r.backbone_id
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

    st.write(f"Backbone: **{_bb_label(rrow)}**")
    st.write(f"Created: {rrow['created_at']}")
    if rrow["note"]:
        st.write(f"Note: {rrow['note']}")
    try:
        params = json.loads(rrow["params_json"] or "{}")
    except Exception:
        params = {}
    if params.get("seq_type"):
        st.caption(f"Sequence type: {params.get('seq_type')}")
    if "remove_stop" in params:
        rs = "Yes" if bool(params.get("remove_stop")) else "No"
        st.caption(f"Remove stop codon (CDS): {rs}")
    if params.get("tm_profile"):
        st.caption(f"Tm profile: {params.get('tm_profile')}")
    if params.get("tm_calc_method"):
        st.caption(f"Tm method: {params.get('tm_calc_method')}")
    if "tm_pair_delta_max" in params:
        st.caption(f"Max |TmF-TmR| (soft): {params.get('tm_pair_delta_max')}")
    seq_type = params.get("seq_type", "")
    remove_stop = bool(params.get("remove_stop")) if "remove_stop" in params else ""

    # Delete controls
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
                conn.execute(
                    f"DELETE FROM primer_runs WHERE user_id=? AND id IN ({placeholders})",
                    [user["id"], *batch_del],
                )
                conn.commit()
                st.success(f"Deleted {len(batch_del)} run(s).")
                st.rerun()

    prof = conn.execute("SELECT * FROM user_profiles WHERE user_id=?", (user["id"],)).fetchone()
    profile = dict(prof) if prof else None

    primers = conn.execute(
        """
        SELECT p.*, g.gene_name
        FROM primers p
        JOIN genes g ON g.id = p.gene_id
        JOIN primer_runs r ON r.id = p.run_id
        WHERE p.run_id=? AND r.user_id=? AND g.user_id=?
        ORDER BY p.id ASC
        """,
        (selected_run, user["id"], user["id"]),
    ).fetchall()

    st.markdown("#### Primers")
    for p in primers[:200]:
        cols = st.columns([3, 2, 2, 2, 2, 3])
        cols[0].write(f"**{p['gene_name']}**")
        cols[1].write(_fmt_tm(p["tm_f"], "TmF"))
        cols[2].write(_fmt_tm(p["tm_r"], "TmR"))
        cols[3].write(f"Amp {p['amplicon_len'] if p['amplicon_len'] is not None else 'NA'} bp")
        cols[4].write(_fmt_delta(p["tm_delta"]))
        cols[5].write(f"Reason: {p['selected_reason'] or 'NA'}")
        st.code(f"{p['primer_name_f']}: {p['final_f']}\n{p['primer_name_r']}: {p['final_r']}", language="text")

    if len(primers) > 200:
        st.info(f"Showing first 200 primer pairs out of {len(primers)} pairs.")

    primer_rows = []
    for p in primers:
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
                "backbone_name": rrow["backbone_name"],
                "gene_name": p["gene_name"],
                "seq_type": seq_type,
                "remove_stop": remove_stop,
                "amplicon_len": p["amplicon_len"],
                "tm_f": p["tm_f"],
                "tm_r": p["tm_r"],
                "tm_delta": p["tm_delta"],
                "tm_profile": p["tm_profile"] or params.get("tm_profile", ""),
                "rescue_tier": p["rescue_tier"] or "",
                "warnings": p["warnings"] or "",
                "selected_reason": p["selected_reason"] or "",
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
            template_bytes = TEMPLATE_PATH.read_bytes()
            template_label = "Built-in template"
        else:
            template_bytes, template_label = get_effective_template_bytes(
                conn, int(user["id"]), int(selected_template_key)
            )
        xbytes = export_from_template_bytes(primer_rows, template_bytes=template_bytes, profile=profile)
        st.download_button(
            "Download Sangon Excel (Template)",
            data=xbytes,
            file_name=f"Sangon_Run_{selected_run}.xlsx",
            key=f"history_download_{selected_run}",
        )
        st.caption(f"Export template: {template_label}")
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
    except Exception as e:
        st.error(f"Re-export blocked by validation: {e}")

    conn.close()
