"""Canonical ordinary-project runner for Primer Project Contract v1."""

from __future__ import annotations

import csv
import io
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from excel_export import export_from_template_bytes
from genbank_export import build_genbank_files
from project_contract import (
    CANONICAL_DIRECTORIES,
    CANONICAL_QC_CHECKS,
    CONSTRUCT_PLAN_FIELDS,
    PRIMER_MASTER_FIELDS,
    PRIMER_QC_FIELDS,
    canonical_project_paths,
)
from project_manifest import resolve_project_manifest, validate_project_manifest
from project_service import (
    ProjectServiceError,
    assemble_project_result,
    normalize_canonical_result,
    normalized_canonical_result_sha256,
)
from sequencing_export import export_sequencing_workbook
from template_resolver import resolve_order_template


class ProjectRunnerError(ValueError):
    """Raised before project artifacts are written when the contract blocks a run."""


def _tsv_bytes(rows: list[dict[str, Any]], fields: Iterable[str]) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "") for field in writer.fieldnames})
    return output.getvalue().encode("utf-8")


def build_project_artifacts(
    manifest: dict[str, Any],
    *,
    sequences: Mapping[str, str],
    identity_audits: Mapping[str, dict[str, Any]],
    structural_qc: Mapping[str, Mapping[str, Any]],
    sequencing_inputs: Mapping[str, Mapping[str, Any]],
    order_template_path: str | Path | None = None,
    sequencing_template_path: str | Path | None = None,
) -> dict[str, Any]:
    """Serialize/export one canonical service result; write nothing on failure."""
    manifest = validate_project_manifest(manifest)
    try:
        canonical = assemble_project_result(
            manifest,
            sequences=sequences,
            identity_audits=identity_audits,
            structural_qc=structural_qc,
            sequencing_inputs=sequencing_inputs,
        )
    except ProjectServiceError as exc:
        raise ProjectRunnerError(str(exc)) from exc
    normalized_canonical = normalize_canonical_result(canonical)

    primer_rows = [row for result in canonical.constructs for row in result.primer_master_rows]
    qc_rows = [row for result in canonical.constructs for row in result.qc_rows]
    plan_rows = [result.construct_plan_row for result in canonical.constructs]
    map_rows = [result.map_input for result in canonical.constructs]
    order_rows = [result.order_input for result in canonical.constructs]
    sequencing_rows = [result.sequencing_input for result in canonical.constructs]
    order_template = resolve_order_template(template_path=order_template_path)
    order_bytes = export_from_template_bytes(
        order_rows,
        template_bytes=order_template.template_bytes,
        purification="HAP",
        nmoles=1,
        tubes=1,
        delivery_format="液体-100μM",
        blank_profile_fields=True,
    )
    sequencing_bytes, sequencing_report = export_sequencing_workbook(
        sequencing_rows,
        explicit_template_path=sequencing_template_path,
        colony_count=canonical.colony_count,
    )
    maps = dict(build_genbank_files(map_rows))
    return {
        "project_id": manifest["project_id"],
        "manifest": json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True).encode("utf-8") + b"\n",
        "construct_plan": _tsv_bytes(plan_rows, CONSTRUCT_PLAN_FIELDS),
        "primer_master": _tsv_bytes(primer_rows, PRIMER_MASTER_FIELDS),
        "primer_qc": _tsv_bytes(qc_rows, PRIMER_QC_FIELDS),
        "maps": maps,
        "order": order_bytes,
        "sequencing": sequencing_bytes,
        "canonical_result": canonical,
        "normalized_canonical_result": normalized_canonical,
        "normalized_canonical_result_sha256": normalized_canonical_result_sha256(
            normalized_canonical
        ),
        "summary": {
            "project": manifest["project_id"], "constructs": len(plan_rows),
            "primers": len(primer_rows), "maps": len(maps),
            "sequencing_samples": sequencing_report["sample_count"],
            "status": "PASS",
        },
    }


def write_project_artifacts(
    project_root: str | Path,
    artifacts: Mapping[str, Any],
    *,
    overwrite_outputs: bool = False,
) -> dict[str, Path]:
    """Reuse/create the control manifest, then write protected derived artifacts."""
    root = Path(project_root)
    if root.name != artifacts["project_id"]:
        raise ProjectRunnerError(
            f"Project directory {root.name!r} must equal project_id {artifacts['project_id']!r}"
        )
    if not (root / "WORK_LOG.md").is_file():
        raise ProjectRunnerError("Project WORK_LOG.md must exist before materializing artifacts")
    paths = canonical_project_paths(root)
    try:
        supplied_manifest = json.loads(bytes(artifacts["manifest"]).decode("utf-8"))
    except (KeyError, TypeError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ProjectRunnerError(f"Artifacts contain an invalid project manifest: {exc}") from exc
    manifest_resolution = resolve_project_manifest(root, supplied_manifest)
    targets = {
        "construct_plan": paths["construct_plan"],
        "primer_master": paths["primer_master"], "primer_qc": paths["primer_qc"],
        "order": paths["order"], "sequencing": paths["sequencing"],
        "summary": paths["report"] / "run_summary.json",
    }
    targets.update({f"map:{name}": paths["maps"] / name for name in artifacts["maps"]})
    existing = [str(path) for path in targets.values() if path.exists()]
    if existing and not overwrite_outputs:
        raise ProjectRunnerError(f"Refusing to overwrite canonical artifacts: {existing[:5]}")
    for relative in CANONICAL_DIRECTORIES:
        (root / relative).mkdir(parents=True, exist_ok=True)
    payloads = {
        "construct_plan": artifacts["construct_plan"],
        "primer_master": artifacts["primer_master"], "primer_qc": artifacts["primer_qc"],
        "order": artifacts["order"], "sequencing": artifacts["sequencing"],
    }
    for key, data in payloads.items():
        targets[key].write_bytes(data)
    for name, data in artifacts["maps"].items():
        targets[f"map:{name}"].write_bytes(data)
    summary_path = targets["summary"]
    summary_path.write_text(
        json.dumps(artifacts["summary"], ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    targets["manifest"] = manifest_resolution.path
    return targets
