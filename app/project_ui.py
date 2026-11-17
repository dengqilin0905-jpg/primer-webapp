"""Streamlit entrypoint for audited Primer Project Contract packages."""

from __future__ import annotations

import streamlit as st
import io
import zipfile
import hashlib

from project_orchestration import build_project_package


def project_contract_page(_user, _cfg) -> None:
    st.subheader("Primer Project Contract")
    st.caption(
        "使用与 Batch/Codex 相同的 manifest、identity preflight 和 shared runner；"
        "生成物先在临时目录完成 conformance audit，再提供 canonical project ZIP。"
    )
    manifest_file = st.file_uploader("project.yaml", type=["yaml", "yml", "json"])
    runtime_file = st.file_uploader("runtime_inputs.json", type=["json"])
    work_log_file = st.file_uploader("WORK_LOG.md", type=["md"])
    state_key='contract_artifact_'+str(_user.get('id','anonymous'))
    inputs=[f.getvalue() if f is not None else None for f in (manifest_file,runtime_file,work_log_file)]
    token=hashlib.sha256(b''.join(len(x or b'').to_bytes(8,'big')+(x or b'') for x in inputs)).hexdigest()
    cached=st.session_state.get(state_key)
    if cached is not None and cached['token']!=token:
        st.session_state.pop(state_key,None)
    if st.button("Build canonical project package"):
        if not all((manifest_file, runtime_file, work_log_file)):
            st.error("请同时提供 project.yaml、runtime_inputs.json 和 WORK_LOG.md。")
            return
        try:
            package, report = build_project_package(
                manifest_file.getvalue(),
                runtime_file.getvalue(),
                work_log_file.getvalue(),
            )
        except Exception as exc:
            st.error(f"项目 contract 校验或生成失败：{exc}")
            return
        st.session_state[state_key]={'token':token,'package':package,'report':report}
    cached=st.session_state.get(state_key)
    if cached is not None and cached['token']==token:
        package,report=cached['package'],cached['report']
        project_id = report["summary"]["project"]
        st.success("Canonical outputs 与 conformance validation 均已通过。")
        st.json(
            {
                "summary": report["summary"],
                "conformance": report["conformance"],
                "normalized_canonical_result_sha256": report[
                    "normalized_canonical_result_sha256"
                ],
            }
        )
        st.download_button(
            "Download canonical project ZIP",
            data=package,
            file_name=f"{project_id}.zip",
            mime="application/zip",
        )

        with zipfile.ZipFile(io.BytesIO(package)) as archive:
            sequencing = archive.read("results/sequencing/sequencing_order.xlsx")
        st.download_button("Download sequencing Excel", sequencing, file_name="sequencing_order.xlsx",
                           mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
