import json
import time

import streamlit as st

from db import can_delete_backbone, get_conn
from file_safety import validate_upload
from primer_design import normalize_seq
from seq_io import parse_uploaded_file

ALLOWED_VECTOR_EXTS = {".fasta", ".fa", ".fna", ".txt", ".gb", ".gbk", ".genbank", ".dna"}


def _snapshot(
    name: str,
    enzyme_sites: str,
    homo_f: str,
    homo_r: str,
    vector_sequence: str,
    insertion_mode: str,
    restriction_site_f: str,
    restriction_site_r: str,
    junction_extra_bases_f: int,
    junction_extra_bases_r: int,
    vector_map_note: str,
) -> dict:
    return {
        "name": (name or "").strip(),
        "enzyme_sites": (enzyme_sites or "").strip(),
        "homo_f": normalize_seq(homo_f or ""),
        "homo_r": normalize_seq(homo_r or ""),
        "vector_sequence": normalize_seq(vector_sequence or ""),
        "insertion_mode": (insertion_mode or "linear_join").strip(),
        "restriction_site_f": normalize_seq(restriction_site_f or ""),
        "restriction_site_r": normalize_seq(restriction_site_r or ""),
        "junction_extra_bases_f": int(junction_extra_bases_f or 0),
        "junction_extra_bases_r": int(junction_extra_bases_r or 0),
        "vector_map_note": (vector_map_note or "").strip(),
    }


def _bb_label(row) -> str:
    enz = (row["enzyme_sites"] or "").strip()
    return f"{row['name']} [{enz}]" if enz else f"{row['name']} [no-enzyme]"


def _load_vector_from_upload(uploaded_file, max_upload_size_mb: int) -> tuple[str, str]:
    validate_upload(uploaded_file, ALLOWED_VECTOR_EXTS, int(max_upload_size_mb))
    recs = parse_uploaded_file(uploaded_file)
    if not recs:
        raise ValueError("No sequence found in uploaded vector file.")
    seq = normalize_seq(recs[0][1])
    if not seq:
        raise ValueError("Uploaded vector sequence is empty after normalization.")
    if len(recs) > 1:
        return seq, f"Loaded first record from file ({len(recs)} records found)."
    return seq, "Vector sequence loaded from uploaded file."


def backbones_page(user, config):
    st.subheader("Backbones")
    conn = get_conn()
    rows = conn.execute(
        "SELECT * FROM backbones WHERE user_id=? ORDER BY name, enzyme_sites, id",
        (user["id"],),
    ).fetchall()

    st.markdown("#### Add / Update backbone")
    name = st.text_input("Backbone name", placeholder="e.g., pGreenII", key="bb_name")
    enzyme_sites = st.text_input("Enzyme sites label (optional)", placeholder="e.g., BamHI/EcoRI", key="bb_enzyme_sites")
    homo_f = st.text_area("HOMO_F (5'->3')", height=80, key="bb_homo_f")
    homo_r = st.text_area("HOMO_R (5'->3')", height=80, key="bb_homo_r")

    st.markdown("##### Vector map settings (optional)")
    vector_upload = st.file_uploader(
        "Upload vector sequence file (.fasta/.gbk/.dna etc.)",
        type=["fasta", "fa", "fna", "txt", "gb", "gbk", "genbank", "dna"],
        accept_multiple_files=False,
        key="bb_vector_upload",
    )
    if st.button("Load vector from uploaded file", key="bb_vector_load_btn"):
        if not vector_upload:
            st.warning("Please choose a vector file first.")
        else:
            try:
                seq, msg = _load_vector_from_upload(vector_upload, int(config.max_upload_size_mb))
                st.session_state["bb_vector_sequence"] = seq
                st.success(f"{msg} Length={len(seq)} bp.")
            except Exception as exc:
                st.error(f"Vector import failed: {exc}")

    vector_sequence = st.text_area(
        "vector_sequence (FASTA/raw DNA accepted)",
        height=120,
        key="bb_vector_sequence",
        help="Paste full vector sequence or load from file above.",
    )
    insertion_mode = st.selectbox(
        "insertion_mode",
        ["linear_join", "restriction_site_join"],
        key="bb_insertion_mode",
    )
    restriction_site_f = st.text_input("restriction_site_F (sequence, optional)", key="bb_restriction_site_f")
    restriction_site_r = st.text_input("restriction_site_R (sequence, optional)", key="bb_restriction_site_r")
    junction_extra_bases_f = st.number_input("junction extra bases (5')", min_value=0, step=1, value=0, key="bb_jextra_f")
    junction_extra_bases_r = st.number_input("junction extra bases (3')", min_value=0, step=1, value=0, key="bb_jextra_r")
    vector_map_note = st.text_area("vector_map_note", height=60, key="bb_vector_map_note")

    curr_snapshot = _snapshot(
        name=name,
        enzyme_sites=enzyme_sites,
        homo_f=homo_f,
        homo_r=homo_r,
        vector_sequence=vector_sequence,
        insertion_mode=insertion_mode,
        restriction_site_f=restriction_site_f,
        restriction_site_r=restriction_site_r,
        junction_extra_bases_f=int(junction_extra_bases_f),
        junction_extra_bases_r=int(junction_extra_bases_r),
        vector_map_note=vector_map_note,
    )
    saved_snapshot = st.session_state.get("bb_saved_snapshot", {})
    if curr_snapshot != saved_snapshot and any(curr_snapshot.values()):
        st.warning("You have unsaved backbone edits.")

    missing = []
    if not curr_snapshot["name"]:
        missing.append("Backbone name")
    if not curr_snapshot["homo_f"]:
        missing.append("HOMO_F")
    if not curr_snapshot["homo_r"]:
        missing.append("HOMO_R")
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
            st.error("Backbone name + HOMO_F + HOMO_R are required.")
        else:
            try:
                conn.execute(
                    """
                    INSERT INTO backbones(
                        user_id, name, enzyme_sites, homo_f, homo_r,
                        vector_sequence, insertion_mode, restriction_site_f, restriction_site_r,
                        junction_extra_bases_f, junction_extra_bases_r, vector_map_note,
                        defaults_json, created_at
                    )
                    VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                    ON CONFLICT(user_id, name, enzyme_sites) DO UPDATE SET
                        homo_f=excluded.homo_f,
                        homo_r=excluded.homo_r,
                        vector_sequence=excluded.vector_sequence,
                        insertion_mode=excluded.insertion_mode,
                        restriction_site_f=excluded.restriction_site_f,
                        restriction_site_r=excluded.restriction_site_r,
                        junction_extra_bases_f=excluded.junction_extra_bases_f,
                        junction_extra_bases_r=excluded.junction_extra_bases_r,
                        vector_map_note=excluded.vector_map_note,
                        defaults_json=excluded.defaults_json
                    """,
                    (
                        user["id"],
                        curr_snapshot["name"],
                        curr_snapshot["enzyme_sites"],
                        curr_snapshot["homo_f"],
                        curr_snapshot["homo_r"],
                        curr_snapshot["vector_sequence"],
                        curr_snapshot["insertion_mode"],
                        curr_snapshot["restriction_site_f"],
                        curr_snapshot["restriction_site_r"],
                        curr_snapshot["junction_extra_bases_f"],
                        curr_snapshot["junction_extra_bases_r"],
                        curr_snapshot["vector_map_note"],
                        json.dumps(defaults, ensure_ascii=False),
                        time.strftime("%Y-%m-%d %H:%M:%S"),
                    ),
                )
                conn.commit()
                st.session_state["bb_saved_snapshot"] = curr_snapshot
                st.success(f"Saved backbone: {curr_snapshot['name']} [{curr_snapshot['enzyme_sites'] or 'no-enzyme'}]")
                st.rerun()
            except Exception as exc:
                st.error(f"Save failed: {exc}")

    st.markdown("#### Your backbones")
    if rows:
        show_rows = []
        for r in rows:
            x = dict(r)
            x["display"] = _bb_label(r)
            x["vector_len"] = len(normalize_seq(r["vector_sequence"] or ""))
            show_rows.append(
                {
                    "id": x["id"],
                    "display": x["display"],
                    "insertion_mode": x.get("insertion_mode", "linear_join"),
                    "restriction_site_f": x.get("restriction_site_f", ""),
                    "restriction_site_r": x.get("restriction_site_r", ""),
                    "vector_len": x["vector_len"],
                    "vector_map_note": x.get("vector_map_note", ""),
                }
            )
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
