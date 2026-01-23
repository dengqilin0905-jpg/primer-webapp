import time, json, re
import streamlit as st
from db import get_conn
from primer_design import design_primers_for_gene
from excel_export import export_sangon_template_bytes

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

def design_page(user):
    st.subheader("Design (Mapping-based batch)")
    conn = get_conn()

    prof = conn.execute("SELECT * FROM user_profiles WHERE user_id=?", (user["id"],)).fetchone()
    profile = dict(prof) if prof else None

    backbones = conn.execute("SELECT * FROM backbones WHERE user_id=? ORDER BY name", (user["id"],)).fetchall()
    if not backbones:
        st.warning("Please add at least one backbone first.")
        conn.close()
        return

    defaults = json.loads((backbones[0]["defaults_json"] or "{}"))

    st.markdown("#### Parameters")
    col1, col2, col3 = st.columns(3)
    with col1:
        remove_stop = st.checkbox("Remove terminal stop codon", value=bool(defaults.get("remove_stop", True)))
        min_len = st.number_input("Min primer length", value=int(defaults.get("min_len", 18)), step=1)
        max_len = st.number_input("Max primer length", value=int(defaults.get("max_len", 30)), step=1)
    with col2:
        tm_min = st.number_input("Tm min (°C)", value=float(defaults.get("tm_min", 55.0)), step=0.5)
        tm_max = st.number_input("Tm max (°C)", value=float(defaults.get("tm_max", 65.0)), step=0.5)
        target_tm = st.number_input("Target Tm (°C)", value=float(defaults.get("target_tm", 60.0)), step=0.5)
    with col3:
        Na = st.number_input("Na+ (mM)", value=50.0, step=5.0)
        Mg = st.number_input("Mg2+ (mM)", value=1.5, step=0.5)
        dNTPs = st.number_input("dNTPs (mM)", value=0.2, step=0.1)
        primer_conc = st.number_input("Primer conc (nM)", value=250.0, step=50.0)

    st.markdown("#### Build mapping (Backbone → Genes)")
    search = st.text_input("Filter genes by name (optional)")
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
        st.session_state["mapping_cart"] = []  # list of {backbone_name, gene_names}

    bname = st.selectbox("Select backbone (one at a time)", [b["name"] for b in backbones])
    selected_genes = st.multiselect("Select genes for this backbone", [g["gene_name"] for g in genes])

    c1, c2, c3 = st.columns([1,1,2])
    with c1:
        if st.button("Add to batch"):
            if not selected_genes:
                st.warning("Please select at least one gene.")
            else:
                st.session_state["mapping_cart"].append({"backbone_name": bname, "gene_names": selected_genes})
                st.success(f"Added: {bname} → {len(selected_genes)} genes")
    with c2:
        if st.button("Clear batch"):
            st.session_state["mapping_cart"] = []
            st.success("Cleared.")
    with c3:
        st.caption("You fully control pairing. No forced Cartesian product.")

    if st.session_state["mapping_cart"]:
        st.markdown("##### Current batch mapping")
        for idx, item in enumerate(st.session_state["mapping_cart"], start=1):
            cols = st.columns([3,6,1])
            cols[0].write(f"{idx}. **{item['backbone_name']}**")
            cols[1].write(", ".join(item["gene_names"]))
            if cols[2].button("X", key=f"rm_{idx}"):
                st.session_state["mapping_cart"].pop(idx-1)
                st.rerun()
    else:
        st.info("Batch is empty. Add at least one mapping above.")

    note = st.text_input("Run note (optional)")

    if st.button("Generate primers for batch mapping"):
        if not st.session_state["mapping_cart"]:
            st.warning("Batch is empty.")
        else:
            group_tag = time.strftime("%Y%m%d_%H%M%S")
            base_note = (note or "").strip()
            run_ids = []
            all_results = []

            params = dict(
                remove_stop=bool(remove_stop),
                min_len=int(min_len),
                max_len=int(max_len),
                tm_min=float(tm_min),
                tm_max=float(tm_max),
                target_tm=float(target_tm),
                Na=float(Na),
                Mg=float(Mg),
                dNTPs=float(dNTPs),
                primer_conc=float(primer_conc),
            )

            gene_by_name = {g["gene_name"]: g for g in genes}

            total_steps = sum(len(it["gene_names"]) for it in st.session_state["mapping_cart"])
            total_steps = max(1, total_steps)
            done = 0
            prog = st.progress(0.0)

            for item in st.session_state["mapping_cart"]:
                bb = next(b for b in backbones if b["name"] == item["backbone_name"])
                gene_list = [gene_by_name[n] for n in item["gene_names"] if n in gene_by_name]

                cur = conn.cursor()
                run_note = f"[group:{group_tag}] " + (base_note if base_note else "")
                cur.execute(
                    "INSERT INTO primer_runs(user_id, backbone_id, params_json, created_at, note) VALUES(?,?,?,?,?)",
                    (user["id"], bb["id"], json.dumps(params, ensure_ascii=False), time.strftime("%Y-%m-%d %H:%M:%S"), run_note or None),
                )
                run_id = cur.lastrowid
                conn.commit()
                run_ids.append(run_id)

                for g in gene_list:
                    out = design_primers_for_gene(
                        gene_name=g["gene_name"],
                        seq=g["sequence"],
                        homo_f=bb["homo_f"],
                        homo_r=bb["homo_r"],
                        remove_stop=bool(remove_stop),
                        min_len=int(min_len),
                        max_len=int(max_len),
                        tm_min=float(tm_min),
                        tm_max=float(tm_max),
                        target_tm=float(target_tm),
                        Na=float(Na),
                        Mg=float(Mg),
                        dNTPs=float(dNTPs),
                        primer_conc=float(primer_conc),
                    )
                    out["primer_name_f"] = _primer_name(bb["name"], g["gene_name"], "F")
                    out["primer_name_r"] = _primer_name(bb["name"], g["gene_name"], "R")
                    out["gene_id"] = g["id"]
                    all_results.append(out)

                    conn.execute(
                        '''
                        INSERT INTO primers(run_id, gene_id, primer_name_f, primer_name_r,
                                            gene_bind_f, gene_bind_r, final_f, final_r,
                                            tm_f, tm_r, gc_f, gc_r, warnings)
                        VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)
                        ''',
                        (
                            run_id, g["id"], out["primer_name_f"], out["primer_name_r"],
                            out["gene_bind_f"], out["gene_bind_r"], out["final_f"], out["final_r"],
                            out.get("tm_f"), out.get("tm_r"), out.get("gc_f"), out.get("gc_r"),
                            out.get("warnings",""),
                        ),
                    )
                    conn.commit()
                    done += 1
                    prog.progress(done / total_steps)

            st.success(f"Generated primers for your mapping. Saved runs: {run_ids}")
            st.session_state["last_results"] = all_results

    if "last_results" in st.session_state and st.session_state["last_results"]:
        st.divider()
        st.markdown("### Latest results (this session)")
        st.dataframe(
            [
                {k: v for k, v in r.items() if k in ["primer_name_f","primer_name_r","tm_f","tm_r","gc_f","gc_r","warnings"]}
                for r in st.session_state["last_results"]
            ],
            use_container_width=True,
        )
        xbytes = export_sangon_template_bytes(st.session_state["last_results"], profile=profile)
        st.download_button("Download Sangon Excel (Template)", data=xbytes, file_name="Sangon_Primers_Template.xlsx")

    conn.close()
