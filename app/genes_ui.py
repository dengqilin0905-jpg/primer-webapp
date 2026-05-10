import time
import streamlit as st
from db import get_conn, can_delete_gene_ids
from seq_io import parse_uploaded_file, parse_fasta_text, derive_gene_name, file_base_name
from file_safety import validate_upload


ALLOWED_GENE_EXTS = {".fasta", ".fa", ".fna", ".txt", ".gb", ".gbk", ".genbank", ".dna"}


def _unique_name(base_name: str, reserved: set[str]) -> str:
    base = (base_name or "").strip() or "Gene"
    if base not in reserved:
        reserved.add(base)
        return base
    i = 2
    while True:
        candidate = f"{base}_{i}"
        if candidate not in reserved:
            reserved.add(candidate)
            return candidate
        i += 1


def genes_page(user, config):
    st.subheader("Genes (Batch import)")
    st.caption("Supported: FASTA/TXT, GenBank (.gb/.gbk/.genbank), SnapGene (.dna).")

    conn = get_conn()

    with st.expander("Utilities", expanded=False):
        if st.button("Normalize existing gene names (prefixTF#)", key="genes_normalize_names"):
            rows = conn.execute("SELECT id, gene_name FROM genes WHERE user_id=?", (user["id"],)).fetchall()
            changed = 0
            for r in rows:
                new_name = derive_gene_name(r["gene_name"])
                if new_name and new_name != r["gene_name"]:
                    conn.execute("UPDATE genes SET gene_name=? WHERE id=?", (new_name, r["id"]))
                    changed += 1
            conn.commit()
            st.success(f"Renamed {changed} gene(s).")

    st.markdown("#### Upload files")
    files = st.file_uploader(
        "Upload sequence files",
        type=["fasta", "fa", "fna", "txt", "gb", "gbk", "genbank", "dna"],
        accept_multiple_files=True,
        key="genes_upload_files",
    )

    st.markdown("#### Paste FASTA / sequence (optional)")
    pasted = st.text_area("Paste FASTA or raw sequence here", height=140, key="genes_pasted_text")

    if st.button("Import", key="genes_import_button"):
        imported = 0
        skipped = 0
        failed_files = []
        import_mapping = []

        existing_names = {
            r["gene_name"]
            for r in conn.execute("SELECT gene_name FROM genes WHERE user_id=?", (user["id"],)).fetchall()
        }

        if pasted.strip():
            recs = parse_fasta_text(pasted)
            for gene_name, seq in recs:
                gene_name = _unique_name(derive_gene_name(gene_name), existing_names)
                if not gene_name or not seq:
                    continue
                try:
                    conn.execute(
                        "INSERT INTO genes(user_id, gene_name, sequence, created_at) VALUES(?,?,?,?)",
                        (user["id"], gene_name, seq, time.strftime("%Y-%m-%d %H:%M:%S")),
                    )
                    imported += 1
                    import_mapping.append({"source": "Pasted input", "gene_name": gene_name, "bp": len(seq)})
                except Exception:
                    skipped += 1
            conn.commit()

        if files:
            for f in files:
                try:
                    validate_upload(f, ALLOWED_GENE_EXTS, int(config.max_upload_size_mb))
                    recs = parse_uploaded_file(f)
                except Exception as e:
                    failed_files.append((f.name, str(e)))
                    continue

                base = derive_gene_name(file_base_name(f.name))
                for idx, (_, seq) in enumerate(recs, start=1):
                    candidate = base if len(recs) == 1 else f"{base}_{idx}"
                    gene_name = _unique_name(candidate, existing_names)
                    if not gene_name or not seq:
                        continue
                    try:
                        conn.execute(
                            "INSERT INTO genes(user_id, gene_name, sequence, created_at) VALUES(?,?,?,?)",
                            (user["id"], gene_name, seq, time.strftime("%Y-%m-%d %H:%M:%S")),
                        )
                        imported += 1
                        import_mapping.append({"source": f.name, "gene_name": gene_name, "bp": len(seq)})
                    except Exception:
                        skipped += 1
                conn.commit()

        st.success(f"Imported {imported} records (skipped duplicates: {skipped}).")
        if import_mapping:
            st.markdown("##### Import mapping (source -> gene name -> bp)")
            st.dataframe(import_mapping, use_container_width=True)
        if failed_files:
            st.warning("Some files could not be parsed.")
            for name, msg in failed_files:
                st.write(f"- {name}: {msg}")

    st.markdown("#### Manage genes")
    gfilter = st.text_input("Filter genes", key="genes_filter")
    if gfilter:
        rows = conn.execute(
            """
            SELECT id, gene_name, LENGTH(sequence) AS len, created_at
            FROM genes
            WHERE user_id=? AND gene_name LIKE ?
            ORDER BY id DESC
            LIMIT 2000
            """,
            (user["id"], f"%{gfilter}%"),
        ).fetchall()
    else:
        rows = conn.execute(
            """
            SELECT id, gene_name, LENGTH(sequence) AS len, created_at
            FROM genes
            WHERE user_id=?
            ORDER BY id DESC
            LIMIT 2000
            """,
            (user["id"],),
        ).fetchall()

    if rows:
        data = [dict(r) for r in rows]
        st.dataframe(data, use_container_width=True)

        gene_options = [r["gene_name"] for r in rows]
        selected_gene_names = st.multiselect(
            "Select genes to delete",
            gene_options,
            key="genes_delete_select",
        )
        confirm_delete = st.checkbox("I understand deletion cannot be undone.", key="genes_delete_confirm")
        if st.button("Delete selected genes", key="genes_delete_button"):
            if not selected_gene_names:
                st.warning("No genes selected.")
            elif not confirm_delete:
                st.warning("Please confirm deletion first.")
            else:
                selected_rows = [r for r in rows if r["gene_name"] in selected_gene_names]
                selected_ids = [int(r["id"]) for r in selected_rows]
                can_delete, reason, refs = can_delete_gene_ids(conn, int(user["id"]), selected_ids)
                if not can_delete:
                    st.error(reason)
                    id_to_name = {int(r["id"]): r["gene_name"] for r in selected_rows}
                    blocker_lines = []
                    for gid, run_ids in refs.items():
                        blocker_lines.append(f"- {id_to_name.get(gid, gid)} blocked by runs: {', '.join(map(str, run_ids[:10]))}")
                    for line in blocker_lines[:20]:
                        st.write(line)
                else:
                    placeholders = ",".join("?" for _ in selected_ids)
                    conn.execute(
                        f"DELETE FROM genes WHERE user_id=? AND id IN ({placeholders})",
                        [user["id"], *selected_ids],
                    )
                    conn.commit()
                    st.success(f"Deleted {len(selected_ids)} gene(s).")
                    st.rerun()
    else:
        st.info("No genes found.")

    conn.close()
