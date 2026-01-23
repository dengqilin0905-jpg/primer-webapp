import streamlit as st
from db import get_conn
from excel_export import export_sangon_template_bytes

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

def history_page(user):
    st.subheader("History")
    conn = get_conn()

    runs = conn.execute(
        '''
        SELECT r.id, r.created_at, r.note, b.name AS backbone_name
        FROM primer_runs r
        JOIN backbones b ON b.id = r.backbone_id
        WHERE r.user_id=?
        ORDER BY r.id DESC
        LIMIT 50
        ''',
        (user["id"],),
    ).fetchall()

    if not runs:
        st.info("No runs yet.")
        conn.close()
        return

    run_ids = [r["id"] for r in runs]
    selected_run = st.selectbox("Select a run", run_ids, format_func=lambda rid: f"Run {rid}")
    rrow = next(r for r in runs if r["id"] == selected_run)

    st.write(f"Backbone: **{rrow['backbone_name']}**")
    st.write(f"Created: {rrow['created_at']}")
    if rrow["note"]:
        st.write(f"Note: {rrow['note']}")

    prof = conn.execute("SELECT * FROM user_profiles WHERE user_id=?", (user["id"],)).fetchone()
    profile = dict(prof) if prof else None

    primers = conn.execute(
        '''
        SELECT p.*, g.gene_name
        FROM primers p
        JOIN genes g ON g.id = p.gene_id
        WHERE p.run_id=?
        ORDER BY p.id ASC
        ''',
        (selected_run,),
    ).fetchall()

    st.markdown("#### Primers")
    for p in primers[:200]:
        cols = st.columns([3,3,3])
        cols[0].write(f"**{p['gene_name']}**")
        cols[1].write(_fmt_tm(p["tm_f"], "TmF"))
        cols[2].write(_fmt_tm(p["tm_r"], "TmR"))
        st.code(f"{p['primer_name_f']}: {p['final_f']}\n{p['primer_name_r']}: {p['final_r']}", language="text")

    if len(primers) > 200:
        st.info(f"Showing first 200 primer pairs out of {len(primers)} pairs.")

    primer_rows = []
    for p in primers:
        primer_rows.append({
            "primer_name_f": p["primer_name_f"],
            "primer_name_r": p["primer_name_r"],
            "final_f": p["final_f"],
            "final_r": p["final_r"],
            "warnings": p["warnings"] or "",
        })
    xbytes = export_sangon_template_bytes(primer_rows, profile=profile)
    st.download_button("Download Sangon Excel (Template)", data=xbytes, file_name=f"Sangon_Run_{selected_run}.xlsx")

    conn.close()
