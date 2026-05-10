import time
import json
import streamlit as st
from db import get_conn, can_delete_backbone


def _snapshot(name: str, enzyme_sites: str, homo_f: str, homo_r: str) -> dict:
    return {
        "name": (name or "").strip(),
        "enzyme_sites": (enzyme_sites or "").strip(),
        "homo_f": (homo_f or "").strip(),
        "homo_r": (homo_r or "").strip(),
    }


def _bb_label(row) -> str:
    enz = (row["enzyme_sites"] or "").strip()
    return f"{row['name']} [{enz}]" if enz else f"{row['name']} [no-enzyme]"


def backbones_page(user):
    st.subheader("Backbones")

    conn = get_conn()
    rows = conn.execute(
        "SELECT * FROM backbones WHERE user_id=? ORDER BY name, enzyme_sites, id",
        (user["id"],),
    ).fetchall()

    st.markdown("#### Add / Update backbone")
    name = st.text_input("Backbone name", placeholder="e.g., pGreenII", key="bb_name")
    enzyme_sites = st.text_input("Enzyme sites (optional)", placeholder="e.g., BamHI/EcoRI", key="bb_enzyme_sites")
    homo_f = st.text_area("Homology arm F (5'->3')", height=90, key="bb_homo_f")
    homo_r = st.text_area("Homology arm R (5'->3')", height=90, key="bb_homo_r")

    curr_snapshot = _snapshot(name, enzyme_sites, homo_f, homo_r)
    saved_snapshot = st.session_state.get("bb_saved_snapshot", {})
    if curr_snapshot != saved_snapshot and any(curr_snapshot.values()):
        st.warning("You have unsaved backbone edits.")

    missing = []
    if not curr_snapshot["name"]:
        missing.append("Backbone name")
    if not curr_snapshot["homo_f"]:
        missing.append("Homology arm F")
    if not curr_snapshot["homo_r"]:
        missing.append("Homology arm R")
    if missing:
        st.caption("Required: " + ", ".join(missing))

    defaults = {
        "remove_stop": True,
        "min_len": 18,
        "max_len": 30,
        "tm_min": 55.0,
        "tm_max": 65.0,
        "target_tm": 60.0,
    }
    if st.button("Save backbone", key="bb_save_button"):
        if missing:
            st.error("Backbone name + homo_f + homo_r are required.")
        else:
            try:
                conn.execute(
                    """
                    INSERT INTO backbones(user_id, name, enzyme_sites, homo_f, homo_r, defaults_json, created_at)
                    VALUES(?,?,?,?,?,?,?)
                    ON CONFLICT(user_id, name, enzyme_sites) DO UPDATE SET
                        homo_f=excluded.homo_f,
                        homo_r=excluded.homo_r,
                        defaults_json=excluded.defaults_json
                    """,
                    (
                        user["id"],
                        curr_snapshot["name"],
                        curr_snapshot["enzyme_sites"],
                        curr_snapshot["homo_f"],
                        curr_snapshot["homo_r"],
                        json.dumps(defaults, ensure_ascii=False),
                        time.strftime("%Y-%m-%d %H:%M:%S"),
                    ),
                )
                conn.commit()
                st.session_state["bb_saved_snapshot"] = curr_snapshot
                st.success(f"Saved backbone: {curr_snapshot['name']} [{curr_snapshot['enzyme_sites'] or 'no-enzyme'}]")
                st.rerun()
            except Exception as e:
                st.error(f"Save failed: {e}")

    st.markdown("#### Your backbones")
    if rows:
        show_rows = []
        for r in rows:
            x = dict(r)
            x["display"] = _bb_label(r)
            show_rows.append(x)
        st.dataframe(show_rows, use_container_width=True)

        labels = [_bb_label(r) for r in rows]
        selected_label = st.selectbox("Select a backbone to delete", labels, key="bb_delete_pick")
        confirm_delete = st.checkbox("I understand deletion cannot be undone.", key="bb_delete_confirm")
        if st.button("Delete selected backbone", key="bb_delete_button"):
            if not confirm_delete:
                st.warning("Please confirm deletion first.")
            else:
                idx = labels.index(selected_label)
                target = rows[idx]
                can_delete, reason, run_ids = can_delete_backbone(conn, int(user["id"]), int(target["id"]))
                if not can_delete:
                    st.error(reason)
                    st.write(f"- {selected_label} blocked by runs: {', '.join(map(str, run_ids[:20]))}")
                else:
                    conn.execute("DELETE FROM backbones WHERE user_id=? AND id=?", (user["id"], target["id"]))
                    conn.commit()
                    st.success(f"Deleted backbone: {selected_label}")
                    st.rerun()
    else:
        st.info("No backbones yet.")
    conn.close()
