import time
import streamlit as st
from db import get_conn
from seq_io import parse_uploaded_file, parse_fasta_text, derive_gene_name

def genes_page(user):
    st.subheader("Genes (Batch import)")
    st.caption("Supported: FASTA/TXT, GenBank (.gb/.gbk). SnapGene .dna is best-effort; recommended workflow is to export FASTA from SnapGene.")

    conn = get_conn()

    with st.expander("Utilities", expanded=False):
        if st.button("Normalize existing gene names (prefixTF#)"):
            rows = conn.execute("SELECT id, gene_name FROM genes WHERE user_id=?", (user['id'],)).fetchall()
            changed = 0
            for r in rows:
                new_name = derive_gene_name(r['gene_name'])
                if new_name and new_name != r['gene_name']:
                    conn.execute("UPDATE genes SET gene_name=? WHERE id=?", (new_name, r['id']))
                    changed += 1
            conn.commit()
            st.success(f"Renamed {changed} gene(s).")

    st.markdown("#### Upload files")
    files = st.file_uploader("Upload sequence files", type=["fasta","fa","fna","txt","gb","gbk","dna"], accept_multiple_files=True)

    st.markdown("#### Paste FASTA / sequence (optional)")
    pasted = st.text_area("Paste FASTA or raw sequence here", height=140)

    if st.button("Import"):
        imported = 0
        skipped = 0
        failed_files = []

        # 1) pasted
        if pasted.strip():
            recs = parse_fasta_text(pasted)
            for gene_name, seq in recs:
                gene_name = gene_name.strip()
                if not gene_name or not seq:
                    continue
                try:
                    conn.execute(
                        "INSERT INTO genes(user_id, gene_name, sequence, created_at) VALUES(?,?,?,?)",
                        (user["id"], gene_name, seq, time.strftime("%Y-%m-%d %H:%M:%S")),
                    )
                    imported += 1
                except Exception:
                    skipped += 1
            conn.commit()

        # 2) files
        if files:
            for f in files:
                try:
                    recs = parse_uploaded_file(f)
                except Exception as e:
                    failed_files.append((f.name, str(e)))
                    continue

                for gene_name, seq in recs:
                    gene_name = gene_name.strip()
                    if not gene_name or not seq:
                        continue
                    try:
                        conn.execute(
                            "INSERT INTO genes(user_id, gene_name, sequence, created_at) VALUES(?,?,?,?)",
                            (user["id"], gene_name, seq, time.strftime("%Y-%m-%d %H:%M:%S")),
                        )
                        imported += 1
                    except Exception:
                        skipped += 1
                conn.commit()

        st.success(f"Imported {imported} records (skipped duplicates: {skipped}).")
        if failed_files:
            st.warning("Some files could not be parsed. Recommended: export FASTA from SnapGene and import FASTA.")
            for name, msg in failed_files:
                st.write(f"- {name}: {msg}")

    st.markdown("#### Your genes (latest 200)")
    rows = conn.execute(
        "SELECT gene_name, LENGTH(sequence) AS len, created_at FROM genes WHERE user_id=? ORDER BY id DESC LIMIT 200",
        (user["id"],),
    ).fetchall()
    if rows:
        st.dataframe([dict(r) for r in rows], use_container_width=True)
    else:
        st.info("No genes yet.")
    conn.close()
