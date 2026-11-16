"""Shared Web and Batch entrypoints for Primer Project Contract projects."""

from __future__ import annotations

import io
import json
import tempfile
import zipfile
from pathlib import Path
from typing import Any, Mapping

from identity_audit import (
    audit_identity,
    collect_identity_evidence,
    discover_active_project_evidence,
    identity_key,
)
from project_contract import PROJECTS_ROOT, HARD_REQUIRED_IDENTITY_SOURCE_KINDS, IDENTITY_PREFLIGHT_SOURCE_KINDS
from project_manifest import (
    load_project_manifest_bytes,
    resolve_project_manifest,
    validate_project_manifest,
)
from project_runner import build_project_artifacts, write_project_artifacts
from project_validator import audit_project


DEFAULT_RUNTIME_INPUT = Path("inputs/runtime_inputs.json")
RUNTIME_INPUT_FIELDS = {
    "sequences",
    "identity_sources",
    "structural_qc",
    "sequencing_inputs",
    "new_identity_approvals",
}


class ProjectOrchestrationError(ValueError):
    """Raised when an entrypoint cannot safely execute the shared runner."""


def _load_json_bytes(data: bytes, *, source: str) -> dict[str, Any]:
    try:
        value = json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ProjectOrchestrationError(f"Unable to read JSON-compatible {source}: {exc}") from exc
    if not isinstance(value, dict):
        raise ProjectOrchestrationError(f"{source} root must be an object")
    return value


def load_runtime_inputs(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    try:
        data = source.read_bytes()
    except OSError as exc:
        raise ProjectOrchestrationError(f"Unable to read runtime inputs {source}: {exc}") from exc
    return validate_runtime_inputs(_load_json_bytes(data, source=str(source)))


def validate_runtime_inputs(data: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(data, Mapping):
        raise ProjectOrchestrationError("runtime inputs root must be an object")
    unknown = sorted(set(data) - RUNTIME_INPUT_FIELDS)
    if unknown:
        raise ProjectOrchestrationError(f"Unsupported runtime input fields: {unknown}")
    if "identity_audits" in data:
        raise ProjectOrchestrationError(
            "Precomputed identity_audits are not accepted; declare identity_sources"
        )
    required = ("sequences", "identity_sources", "structural_qc", "sequencing_inputs")
    missing = [field for field in required if not isinstance(data.get(field), Mapping)]
    if missing:
        raise ProjectOrchestrationError(f"Missing runtime input objects: {missing}")
    sources = data["identity_sources"]
    unknown_kinds = sorted(set(sources) - set(IDENTITY_PREFLIGHT_SOURCE_KINDS))
    if unknown_kinds:
        raise ProjectOrchestrationError(f"Unsupported identity source kinds: {unknown_kinds}")
    normalized_sources: dict[str, list[str]] = {}
    for kind, paths in sources.items():
        if not isinstance(paths, list) or not paths or not all(isinstance(path, str) and path for path in paths):
            raise ProjectOrchestrationError(f"identity_sources.{kind} must be a non-empty path list")
        normalized_sources[kind] = list(paths)
    missing_hard = [kind for kind in HARD_REQUIRED_IDENTITY_SOURCE_KINDS if kind not in normalized_sources]
    if missing_hard:
        raise ProjectOrchestrationError(f"Missing hard-required identity sources: {missing_hard}")
    normalized = dict(data)
    normalized["identity_sources"] = normalized_sources
    approvals = normalized.get("new_identity_approvals", {})
    if not isinstance(approvals, Mapping):
        raise ProjectOrchestrationError("new_identity_approvals must be an object")
    normalized["new_identity_approvals"] = dict(approvals)
    return normalized


def _validate_identity_source_provenance(
    manifest: Mapping[str, Any], sources: Mapping[str, list[str]]
) -> None:
    runtime_paths = sorted(path for paths in sources.values() for path in paths)
    declared_paths = sorted(str(path) for path in manifest["provenance"]["identity_sources"])
    if runtime_paths != declared_paths:
        raise ProjectOrchestrationError(
            "manifest provenance.identity_sources must exactly match runtime identity source paths"
        )


def _approve_new_identities(
    constructs: list[Mapping[str, Any]], audits: dict[str, dict[str, Any]], approvals: Mapping[str, Any], evidence: list[Any]
) -> dict[str, dict[str, Any]]:
    """Approve consecutive new TF identities once, in manifest order."""
    requested_owner: dict[str, str] = {}
    for item in evidence:
        requested_owner.setdefault(identity_key(item.canonical_identity), item.gene_id)
    new_index = 0
    for construct in constructs:
        gene_id = construct["gene_id"]
        audit = audits[gene_id]
        requested = str(construct.get("canonical_identity") or "").strip()
        key = identity_key(requested)
        owner = requested_owner.get(key)
        if owner and owner != gene_id:
            audit.update({"decision": "BLOCK_IDENTITY_COLLISION", "reason": f"{requested} is reserved by {owner}"})
            continue
        if audit["decision"] != "READY_FOR_NEW_IDENTITY":
            continue
        approval = approvals.get(gene_id)
        required = ("canonical_identity", "approved_by", "approved_at")
        if not isinstance(approval, Mapping) or any(not str(approval.get(field) or "").strip() for field in required):
            continue
        approved_identity = str(approval["canonical_identity"]).strip()
        if identity_key(approved_identity) != key:
            raise ProjectOrchestrationError(f"{gene_id}: approved identity conflicts with manifest identity")
        expected = f"TF{int(audit['next_available_tf_id'][2:]) + new_index}"
        if identity_key(approved_identity) != identity_key(expected):
            audit.update({"decision": "BLOCK_IDENTITY_ALLOCATION", "reason": f"Expected consecutive identity {expected}, got {approved_identity}"})
            continue
        if key in requested_owner and requested_owner[key] != gene_id:
            audit.update({"decision": "BLOCK_IDENTITY_COLLISION", "reason": f"Duplicate requested identity {approved_identity}"})
            continue
        audit.update({"decision": "APPROVED_NEW_IDENTITY", "canonical_identity": approved_identity, "approval": dict(approval), "reason": "New identity explicitly approved after master continuity audit"})
        requested_owner[key] = gene_id
        new_index += 1
    return audits


def build_contract_artifacts(
    manifest: Mapping[str, Any],
    runtime_inputs: Mapping[str, Any],
    *,
    order_template_path: str | Path | None = None,
    sequencing_template_path: str | Path | None = None,
    active_projects_root: str | Path | None = None,
    exclude_project_root: str | Path | None = None,
) -> dict[str, Any]:
    """Run source-backed identity preflight, then delegate all outputs to the runner."""
    checked_manifest = validate_project_manifest(dict(manifest))
    checked_inputs = validate_runtime_inputs(runtime_inputs)
    sources = checked_inputs["identity_sources"]
    _validate_identity_source_provenance(checked_manifest, sources)
    evidence = collect_identity_evidence(sources)
    evidence.extend(discover_active_project_evidence(
        projects_root=active_projects_root or PROJECTS_ROOT,
        exclude_project_root=exclude_project_root,
    ))
    approvals = checked_inputs["new_identity_approvals"]
    identity_audits = {}
    for construct in checked_manifest["constructs"]:
        gene_id = construct["gene_id"]
        audit = audit_identity(
            gene_id,
            evidence,
            requested_identity=construct.get("canonical_identity") or None,
        )
        identity_audits[gene_id] = audit
    identity_audits = _approve_new_identities(
        checked_manifest["constructs"], identity_audits, approvals, evidence
    )
    return build_project_artifacts(
        checked_manifest,
        sequences=checked_inputs["sequences"],
        identity_audits=identity_audits,
        structural_qc=checked_inputs["structural_qc"],
        sequencing_inputs=checked_inputs["sequencing_inputs"],
        order_template_path=order_template_path,
        sequencing_template_path=sequencing_template_path,
    )


def run_project_directory(
    project_root: str | Path,
    *,
    manifest: Mapping[str, Any] | None = None,
    runtime_input_path: str | Path | None = None,
    order_template_path: str | Path | None = None,
    sequencing_template_path: str | Path | None = None,
    overwrite_outputs: bool = False,
) -> dict[str, Any]:
    """Batch/Codex entrypoint: load, run, write, and audit one project."""
    root = Path(project_root)
    resolution = resolve_project_manifest(
        root, dict(manifest) if manifest is not None else None
    )
    runtime_path = Path(runtime_input_path) if runtime_input_path else root / DEFAULT_RUNTIME_INPUT
    runtime_inputs = load_runtime_inputs(runtime_path)
    artifacts = build_contract_artifacts(
        resolution.manifest,
        runtime_inputs,
        order_template_path=order_template_path,
        sequencing_template_path=sequencing_template_path,
        active_projects_root=root.parent,
        exclude_project_root=root,
    )
    canonical_runtime = root / DEFAULT_RUNTIME_INPUT
    if canonical_runtime.exists():
        if load_runtime_inputs(canonical_runtime) != runtime_inputs:
            raise ProjectOrchestrationError("Canonical runtime inputs conflict with explicit input")
    elif runtime_path.resolve() != canonical_runtime.resolve():
        canonical_runtime.parent.mkdir(parents=True, exist_ok=True)
        # Preserve the reviewed input bytes; never replace an existing snapshot.
        with canonical_runtime.open("xb") as handle:
            handle.write(runtime_path.read_bytes())
    paths = write_project_artifacts(
        root, artifacts, overwrite_outputs=overwrite_outputs
    )
    report = audit_project(root)
    if report["status"] != "PASS":
        raise ProjectOrchestrationError(
            f"Project outputs failed conformance validation: {report['failed_checks']} checks"
        )
    return {
        "summary": artifacts["summary"],
        "paths": paths,
        "conformance": report,
        "manifest_action": resolution.action,
        "normalized_canonical_result": artifacts["normalized_canonical_result"],
        "normalized_canonical_result_sha256": artifacts[
            "normalized_canonical_result_sha256"
        ],
    }


def build_project_package(
    manifest_bytes: bytes,
    runtime_input_bytes: bytes,
    work_log_bytes: bytes,
    *,
    order_template_path: str | Path | None = None,
    sequencing_template_path: str | Path | None = None,
    active_projects_root: str | Path | None = None,
    exclude_project_root: str | Path | None = None,
) -> tuple[bytes, dict[str, Any]]:
    """Web entrypoint: build an audited canonical project ZIP without server-project writes."""
    manifest = load_project_manifest_bytes(manifest_bytes, source="project.yaml")
    runtime_inputs = validate_runtime_inputs(
        _load_json_bytes(runtime_input_bytes, source="runtime_inputs.json")
    )
    if not work_log_bytes.strip():
        raise ProjectOrchestrationError("WORK_LOG.md must not be empty")
    artifacts = build_contract_artifacts(
        manifest,
        runtime_inputs,
        order_template_path=order_template_path,
        sequencing_template_path=sequencing_template_path,
        active_projects_root=active_projects_root,
        exclude_project_root=exclude_project_root,
    )
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / manifest["project_id"]
        root.mkdir()
        (root / "WORK_LOG.md").write_bytes(work_log_bytes)
        (root / "project.yaml").write_bytes(manifest_bytes)
        write_project_artifacts(root, artifacts)
        (root / DEFAULT_RUNTIME_INPUT).write_bytes(runtime_input_bytes)
        report = audit_project(root)
        if report["status"] != "PASS":
            raise ProjectOrchestrationError(
                f"Web project package failed conformance validation: {report['failed_checks']} checks"
            )
        output = io.BytesIO()
        with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for path in sorted(root.rglob("*")):
                relative = path.relative_to(root)
                if path.is_dir():
                    archive.writestr(f"{relative.as_posix()}/", b"")
                else:
                    archive.write(path, relative.as_posix())
    return output.getvalue(), {
        "summary": artifacts["summary"],
        "conformance": report,
        "normalized_canonical_result": artifacts["normalized_canonical_result"],
        "normalized_canonical_result_sha256": artifacts[
            "normalized_canonical_result_sha256"
        ],
    }
