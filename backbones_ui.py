import time, json
import streamlit as st
from db import get_conn

def backbones_page(user):
    st.subheader("Backbones")

    conn = get_conn()
    rows = conn.execute("SELECT * FROM backbones WHERE user_id=? ORDER BY name", (user["id"],)).fetchall()

    st.markdown("#### Add / Update backbone")
    name = st.text_input("Backbone name", placeholder="e.g., pHW, pCAMBIA1300-FLAG")
    enzyme_sites = st.text_input("Enzyme sites (optional)", placeholder="e.g., BamHI/EcoRI")
    homo_f = st.text_area("Homology arm F (5'->3')", height=90)
    homo_r = st.text_area("Homology arm R (5'->3')", height=90)

    defaults = {
        "remove_stop": True,
        "min_len": 18,
        "max_len": 30,
        "tm_min": 55.0,
        "tm_max": 65.0,
        "target_tm": 60.0,
    }
    if st.button("Save backbone"):
        if not name.strip() or not homo_f.strip() or not homo_r.strip():
            st.error("Backbone name + homo_f + homo_r are required.")
        else:
            try:
                conn.execute(
                    '''
                    INSERT INTO backbones(user_id, name, enzyme_sites, homo_f, homo_r, defaults_json, created_at)
                    VALUES(?,?,?,?,?,?,?)
                    ON CONFLICT(user_id, name) DO UPDATE SET
                        enzyme_sites=excluded.enzyme_sites,
                        homo_f=excluded.homo_f,
                        homo_r=excluded.homo_r,
                        defaults_json=excluded.defaults_json
                    ''',
                    (
                        user["id"],
                        name.strip(),
                        enzyme_sites.strip() or None,
                        homo_f.strip(),
                        homo_r.strip(),
                        json.dumps(defaults, ensure_ascii=False),
                        time.strftime("%Y-%m-%d %H:%M:%S"),
                    ),
                )
                conn.commit()
                st.success("Saved.")
                st.rerun()
            except Exception as e:
                st.error(f"Save failed: {e}")

    st.markdown("#### Your backbones")
    if rows:
        st.dataframe([dict(r) for r in rows], use_container_width=True)
    else:
        st.info("No backbones yet.")
    conn.close()
