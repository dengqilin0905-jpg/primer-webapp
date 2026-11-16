"""Read-only conformance audit for canonical primer projects."""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping

from Bio import SeqIO
from openpyxl import load_workbook

from project_contract import (
    CANONICAL_DIRECTORIES,
    CANONICAL_QC_CHECKS,
    CONSTRUCT_PLAN_FIELDS,
    CURRENT_STRICT,
    LEGACY_COMPATIBILITY,
    PRIMER_MASTER_FIELDS,
    PRIMER_QC_FIELDS,
    REQUIRED_PROJECT_FILES,
    canonical_project_paths,
)
from project_manifest import ProjectManifestError, load_project_manifest
from project_migration import DEFAULT_MIGRATION_REGISTRY, audit_migration_registry, load_migration_registry
from primer_design import normalize_seq, revcomp
from sequencing_export import resolve_sequencing_template
from template_resolver import resolve_order_template


from resource_config import resource_path,state_root
DEFAULT_PRIMER_SKILL_ROOT = resource_path('skill_root',state_root()/'skill')
SEMANTIC_RUNTIME_INPUT = Path("inputs/runtime_inputs.json")


def _report_value(value: Any) -> str:
    if isinstance(value, (dict, list, tuple)):
        rendered = json.dumps(value, ensure_ascii=False, sort_keys=True, default=str)
    else:
        rendered = str(value)
    if len(rendered) > 500:
        return f"{rendered[:500]}... [length={len(rendered)}]"
    return rendered


def _registered_validation_profile(
    root: Path,
    migration_registry: str | Path,
) -> tuple[str, dict[str, Any], str | None]:
    """Resolve compatibility only from controlled registry metadata."""
    try:
        registry = load_migration_registry(migration_registry)
        registered_root = Path(str(registry.get("project_root") or "")).resolve()
        relative = root.resolve().relative_to(registered_root)
    except (Exception, ValueError):
        return CURRENT_STRICT, {}, None
    if len(relative.parts) != 1:
        return CURRENT_STRICT, {}, None
    entry = next(
        (
            item
            for item in registry.get("projects", [])
            if isinstance(item, dict) and item.get("project_id") == relative.name
        ),
        None,
    )
    if not entry or entry.get("validation_profile") != LEGACY_COMPATIBILITY:
        return CURRENT_STRICT, dict(entry or {}), None
    missing = []
    if entry.get("immutable") is not True:
        missing.append("immutable=true")
    for field in ("contract_generation", "annotation_snapshot_provenance"):
        if not str(entry.get(field) or "").strip():
            missing.append(field)
    if missing:
        return CURRENT_STRICT, dict(entry), f"Unsafe legacy registration; missing {missing}"
    return LEGACY_COMPATIBILITY, dict(entry), None


def _read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return tuple(reader.fieldnames or ()), [dict(row) for row in reader]


def _row_strings(rows: Iterable[Mapping[str, Any]], fields: Iterable[str]) -> list[dict[str, str]]:
    return [
        {field: "" if row.get(field) is None else str(row.get(field)) for field in fields}
        for row in rows
    ]


def _compare_table(
    record: Callable[..., None],
    *,
    artifact: str,
    actual_path: Path,
    fields: tuple[str, ...],
    expected_rows: list[dict[str, str]],
    key_fields: tuple[str, ...],
) -> None:
    try:
        actual_fields, actual_rows = _read_tsv(actual_path)
    except Exception as exc:
        record(
            f"semantic:{artifact}:read",
            False,
            str(exc),
            artifact=artifact,
            field="read",
            expected="readable TSV",
            actual=str(exc),
        )
        return
    record(
        f"semantic:{artifact}:schema",
        actual_fields == fields,
        f"columns={list(actual_fields)}",
        artifact=artifact,
        field="columns",
        expected=list(fields),
        actual=list(actual_fields),
    )
    record(
        f"semantic:{artifact}:row_count",
        len(actual_rows) == len(expected_rows),
        f"rows={len(actual_rows)}; expected={len(expected_rows)}",
        artifact=artifact,
        field="row_count",
        expected=len(expected_rows),
        actual=len(actual_rows),
    )

    def key(row: Mapping[str, str]) -> tuple[str, ...]:
        return tuple(str(row.get(field) or "") for field in key_fields)

    actual_index = {key(row): row for row in actual_rows}
    expected_index = {key(row): row for row in expected_rows}
    record(
        f"semantic:{artifact}:primary_keys",
        len(actual_index) == len(actual_rows) and set(actual_index) == set(expected_index),
        f"actual_keys={len(actual_index)}; expected_keys={len(expected_index)}",
        artifact=artifact,
        field="primary_keys",
        expected=sorted(expected_index),
        actual=sorted(actual_index),
    )
    for row_key in sorted(set(actual_index) & set(expected_index)):
        actual = actual_index[row_key]
        expected = expected_index[row_key]
        construct = actual.get("full_construct_name") or actual.get("gene_id") or "/".join(row_key)
        for field in fields:
            if actual.get(field, "") != expected.get(field, ""):
                record(
                    f"semantic:{artifact}:{'/'.join(row_key)}:{field}",
                    False,
                    f"{field} differs",
                    artifact=artifact,
                    construct=construct,
                    field=field,
                    expected=expected.get(field, ""),
                    actual=actual.get(field, ""),
                )


def _parse_genbank(data: bytes):
    return SeqIO.read(io.StringIO(data.decode("utf-8")), "genbank")


def _feature_signature(record, *, strict: bool) -> list[dict[str, Any]]:
    output = []
    for feature in record.features:
        qualifiers = {
            key: [str(value) for value in values]
            for key, values in feature.qualifiers.items()
        }
        if not strict:
            qualifiers = {
                key: qualifiers[key]
                for key in ("label", "gene", "locus_tag", "codon_start")
                if key in qualifiers
            }
        output.append(
            {
                "type": feature.type,
                "location": str(feature.location),
                "qualifiers": qualifiers,
            }
        )
    return output


def _compare_maps(
    record: Callable[..., None],
    *,
    root: Path,
    expected_maps: Mapping[str, bytes],
    strict: bool,
    canonical_result: Any,
) -> None:
    maps_dir = canonical_project_paths(root)["maps"]
    actual_names = sorted(path.name for path in maps_dir.glob("*.gbk")) if maps_dir.is_dir() else []
    expected_names = sorted(expected_maps)
    record(
        "semantic:maps:filenames",
        actual_names == expected_names,
        f"maps={len(actual_names)}; expected={len(expected_names)}",
        artifact="maps",
        field="filenames",
        expected=expected_names,
        actual=actual_names,
    )
    by_name = {
        item.map_input["output_filename"]: item for item in canonical_result.constructs
    }
    for name in sorted(set(actual_names) & set(expected_names)):
        item = by_name[name]
        actual_record = SeqIO.read(maps_dir / name, "genbank")
        expected_record = _parse_genbank(expected_maps[name])
        comparisons = {
            "record_id": (expected_record.id, actual_record.id),
            "record_name": (expected_record.name, actual_record.name),
            "description": (expected_record.description, actual_record.description),
            "sequence_length": (len(expected_record.seq), len(actual_record.seq)),
            "sequence": (str(expected_record.seq).upper(), str(actual_record.seq).upper()),
            "features": (
                _feature_signature(expected_record, strict=strict),
                _feature_signature(actual_record, strict=strict),
            ),
        }
        for field, (expected, actual) in comparisons.items():
            record(
                f"semantic:map:{name}:{field}",
                expected == actual,
                f"map={name}; field={field}",
                artifact=f"results/maps/{name}",
                construct=item.full_construct_name,
                field=field,
                expected=expected,
                actual=actual,
            )
        insert_features = [
            feature
            for feature in actual_record.features
            if (feature.qualifiers.get("locus_tag") or [None])[0] == item.gene_id
        ]
        record(
            f"semantic:map:{name}:insert_identity",
            len(insert_features) == 1,
            f"insert_features={len(insert_features)}",
            artifact=f"results/maps/{name}",
            construct=item.full_construct_name,
            field="insert_locus_tag",
            expected=item.gene_id,
            actual=[feature.qualifiers.get("locus_tag") for feature in insert_features],
        )
        if insert_features:
            feature = insert_features[0]
            insert_sequence = str(feature.extract(actual_record.seq)).upper()
            record(
                f"semantic:map:{name}:insert_sequence",
                insert_sequence == item.sequences.insert,
                f"insert_length={len(insert_sequence)}",
                artifact=f"results/maps/{name}",
                construct=item.full_construct_name,
                field="insert_sequence",
                expected=item.sequences.insert,
                actual=insert_sequence,
            )
            label = (feature.qualifiers.get("label") or [""])[0]
            gene = (feature.qualifiers.get("gene") or [""])[0]
            record(
                f"semantic:map:{name}:canonical_symbol",
                label == item.annotation.resolved_symbol and gene == item.annotation.resolved_symbol,
                f"label={label}; gene={gene}",
                artifact=f"results/maps/{name}",
                construct=item.full_construct_name,
                field="canonical_symbol",
                expected=item.annotation.resolved_symbol,
                actual={"label": label, "gene": gene},
            )
        if strict:
            notes = "\n".join(
                str(value)
                for feature in actual_record.features
                for value in feature.qualifiers.get("note", [])
            )
            for field, expected in (
                ("annotation_source", item.provenance.annotation_source),
                ("annotation_source_path", item.provenance.annotation_source_path),
                ("fallback_type", item.annotation.fallback_type),
            ):
                record(
                    f"semantic:map:{name}:{field}",
                    any(
                        re.sub(r"\s+", "", value) == re.sub(r"\s+", "", str(expected))
                        for value in re.findall(re.escape(field) + r"=([^;]+)", notes)
                    ),
                    f"provenance field {field}",
                    artifact=f"results/maps/{name}",
                    construct=item.full_construct_name,
                    field=field,
                    expected=expected,
                    actual=notes,
                )
        else:
            notes = "\n".join(
                str(value)
                for feature in actual_record.features
                for value in feature.qualifiers.get("note", [])
            )
            for field, expected in (
                ("gene_ID", item.gene_id),
                ("canonical_identity", item.canonical_identity),
            ):
                marker = f"{field}={expected}"
                record(
                    f"semantic:map:{name}:legacy_{field}",
                    marker in notes,
                    marker,
                    artifact=f"results/maps/{name}",
                    construct=item.full_construct_name,
                    field=field,
                    expected=marker,
                    actual=notes,
                )


def _order_rows(data: bytes) -> list[dict[str, Any]]:
    workbook = load_workbook(io.BytesIO(data), data_only=False)
    for sheet in workbook.worksheets:
        for row_number in range(1, min(sheet.max_row, 40) + 1):
            values = [sheet.cell(row_number, column).value for column in range(1, sheet.max_column + 1)]
            labels = [str(value or "").strip() for value in values]
            name_columns = [index + 1 for index, value in enumerate(labels) if "Primer名称" in value]
            sequence_columns = [index + 1 for index, value in enumerate(labels) if "序列" in value and "保密" not in value]
            if not name_columns or not sequence_columns:
                continue
            name_column = name_columns[0]
            sequence_column = sequence_columns[0]
            remark_column = next(
                (index + 1 for index, value in enumerate(labels) if value == "备注"),
                None,
            )
            observed = []
            for row in range(row_number + 1, sheet.max_row + 1):
                name = sheet.cell(row, name_column).value
                sequence = normalize_seq(sheet.cell(row, sequence_column).value or "")
                if name and sequence:
                    observed.append(
                        {
                            "primer_name": str(name).strip(),
                            "sequence": sequence,
                            "primer_length": len(sequence),
                            "remarks": str(sheet.cell(row, remark_column).value or "").strip()
                            if remark_column
                            else "",
                        }
                    )
            return observed
    raise ValueError("Unable to locate primer order header")


def _sequencing_rows(data: bytes) -> list[list[Any]]:
    workbook = load_workbook(io.BytesIO(data), data_only=False)
    sheet = workbook["上传模版"]
    return [
        [sheet.cell(row, column).value for column in range(1, 11)]
        for row in range(4, sheet.max_row + 1)
        if any(sheet.cell(row, column).value is not None for column in range(1, 11))
    ]


def _compare_workbook_rows(
    record: Callable[..., None],
    *,
    artifact: str,
    expected: Any,
    actual: Any,
) -> None:
    record(
        f"semantic:{artifact}:rows",
        expected == actual,
        f"actual_rows={len(actual)}; expected_rows={len(expected)}",
        artifact=artifact,
        field="rows",
        expected=expected,
        actual=actual,
    )


class _LegacySnapshotResolver:
    def __init__(self, manifest: Mapping[str, Any], provenance: str):
        self.symbols = {
            str(item["gene_id"]): str(item["symbol"])
            for item in manifest["constructs"]
        }
        self.provenance = provenance

    def resolve(self, gene_id: str, *, input_hint: str | None = None):
        from rice_annotation import RiceAnnotationResult

        symbol = self.symbols[str(gene_id)]
        return RiceAnnotationResult(
            gene_id=str(gene_id),
            resolved_symbol=symbol,
            official_rice_symbol=None if symbol == "unknown" else symbol,
            at_homolog_id=None,
            at_homolog_symbol=None,
            source="LEGACY_ANNOTATION_SNAPSHOT",
            source_path=self.provenance,
            fallback_type="LEGACY_SNAPSHOT",
            warning="",
            input_hint=input_hint,
            hint_conflict=False,
        )


def _legacy_expected(
    root: Path,
    manifest: Mapping[str, Any],
    registration: Mapping[str, Any],
) -> dict[str, Any]:
    from excel_export import export_from_template_bytes
    from genbank_export import build_genbank_files
    from project_service import assemble_project_result
    from sequencing_export import export_sequencing_workbook

    runtime = json.loads((root / SEMANTIC_RUNTIME_INPUT).read_text(encoding="utf-8"))
    _, plan_rows = _read_tsv(canonical_project_paths(root)["construct_plan"])
    decisions = {row["gene_id"]: row["identity_decision"] for row in plan_rows}
    identity_audits = {
        str(item["gene_id"]): {
            "decision": decisions.get(str(item["gene_id"]), ""),
            "canonical_identity": str(item.get("canonical_identity") or ""),
        }
        for item in manifest["constructs"]
    }
    canonical = assemble_project_result(
        manifest,
        sequences=runtime["sequences"],
        identity_audits=identity_audits,
        structural_qc=runtime["structural_qc"],
        sequencing_inputs=runtime["sequencing_inputs"],
        rice_annotation_resolver=_LegacySnapshotResolver(
            manifest,
            str(registration["annotation_snapshot_provenance"]),
        ),
    )
    plan = [dict(item.construct_plan_row) for item in canonical.constructs]
    primers = [dict(row) for item in canonical.constructs for row in item.primer_master_rows]
    qc = [dict(row) for item in canonical.constructs for row in item.qc_rows]
    for row in qc:
        if row["check"] == "gene_identity":
            symbol = next(
                item.annotation.resolved_symbol
                for item in canonical.constructs
                if item.gene_id == row["gene_id"]
            )
            row["detail"] = f"gene_id={row['gene_id']}; symbol={symbol}"
    map_rows = [dict(item.map_input) for item in canonical.constructs]
    order_rows = [dict(item.order_input) for item in canonical.constructs]
    sequencing_rows = [dict(item.sequencing_input) for item in canonical.constructs]
    order_template = resolve_order_template()
    order = export_from_template_bytes(
        order_rows,
        template_bytes=order_template.template_bytes,
        purification="HAP",
        nmoles=1,
        tubes=1,
        delivery_format="液体-100μM",
        blank_profile_fields=True,
    )
    sequencing, sequencing_report = export_sequencing_workbook(
        sequencing_rows,
        colony_count=canonical.colony_count,
    )
    return {
        "canonical_result": canonical,
        "construct_plan_rows": _row_strings(plan, CONSTRUCT_PLAN_FIELDS),
        "primer_master_rows": _row_strings(primers, PRIMER_MASTER_FIELDS),
        "primer_qc_rows": _row_strings(qc, PRIMER_QC_FIELDS),
        "maps": dict(build_genbank_files(map_rows)),
        "order": order,
        "sequencing": sequencing,
        "summary": {
            "project": manifest["project_id"],
            "constructs": len(plan),
            "primers": len(primers),
            "maps": len(map_rows),
            "sequencing_samples": sequencing_report["sample_count"],
            "status": "PASS",
        },
    }


def _strict_expected(root: Path, manifest: Mapping[str, Any]) -> dict[str, Any]:
    from project_orchestration import build_contract_artifacts, load_runtime_inputs

    runtime = load_runtime_inputs(root / SEMANTIC_RUNTIME_INPUT)
    artifacts = build_contract_artifacts(
        manifest,
        runtime,
        active_projects_root=root.parent,
        exclude_project_root=root,
    )
    for key, fields in (
        ("construct_plan", CONSTRUCT_PLAN_FIELDS),
        ("primer_master", PRIMER_MASTER_FIELDS),
        ("primer_qc", PRIMER_QC_FIELDS),
    ):
        reader = csv.DictReader(io.StringIO(artifacts[key].decode("utf-8")), delimiter="\t")
        artifacts[f"{key}_rows"] = [dict(row) for row in reader]
        artifacts[f"{key}_fields"] = tuple(reader.fieldnames or fields)
    return artifacts


def _audit_semantics(
    root: Path,
    manifest: Mapping[str, Any],
    profile: str,
    registration: Mapping[str, Any],
    record: Callable[..., None],
) -> None:
    try:
        expected = (
            _legacy_expected(root, manifest, registration)
            if profile == LEGACY_COMPATIBILITY
            else _strict_expected(root, manifest)
        )
    except Exception as exc:
        record(
            "semantic:canonical_expected",
            False,
            str(exc),
            artifact="inputs/runtime_inputs.json",
            field="canonical_expected",
            expected="reproducible canonical result",
            actual=str(exc),
        )
        return

    paths = canonical_project_paths(root)
    for artifact, path_key, fields, keys in (
        ("construct_plan", "construct_plan", CONSTRUCT_PLAN_FIELDS, ("gene_id",)),
        (
            "primer_master",
            "primer_master",
            PRIMER_MASTER_FIELDS,
            ("gene_id", "primer_direction"),
        ),
        ("primer_qc", "primer_qc", PRIMER_QC_FIELDS, ("gene_id", "check")),
    ):
        _compare_table(
            record,
            artifact=artifact,
            actual_path=paths[path_key],
            fields=fields,
            expected_rows=expected[f"{artifact}_rows"],
            key_fields=keys,
        )

    canonical = expected["canonical_result"]
    record(
        "semantic:construct_primer_relationship",
        len(canonical.constructs) * 2
        == sum(len(item.primers) for item in canonical.constructs),
        f"constructs={len(canonical.constructs)}; primers={sum(len(item.primers) for item in canonical.constructs)}",
        artifact="primer_master",
        field="construct_primer_count",
        expected=len(canonical.constructs) * 2,
        actual=sum(len(item.primers) for item in canonical.constructs),
    )
    record(
        "semantic:canonical_qc_relationship",
        all(
            {row["check"] for row in item.qc_rows} == set(CANONICAL_QC_CHECKS)
            for item in canonical.constructs
        ),
        f"checks_per_construct={len(CANONICAL_QC_CHECKS)}",
        artifact="primer_qc",
        field="canonical_qc_checks",
        expected=list(CANONICAL_QC_CHECKS),
        actual=[sorted(row["check"] for row in item.qc_rows) for item in canonical.constructs],
    )
    _compare_maps(
        record,
        root=root,
        expected_maps=expected["maps"],
        strict=profile == CURRENT_STRICT,
        canonical_result=canonical,
    )
    try:
        actual_order = _order_rows(paths["order"].read_bytes())
        expected_order = _order_rows(expected["order"])
    except Exception as exc:
        record(
            "semantic:order:read",
            False,
            str(exc),
            artifact="results/order/primer_order.xlsx",
            field="read",
            expected="readable canonical order workbook",
            actual=str(exc),
        )
    else:
        _compare_workbook_rows(
            record,
            artifact="results/order/primer_order.xlsx",
            expected=expected_order,
            actual=actual_order,
        )
    try:
        actual_sequencing = _sequencing_rows(paths["sequencing"].read_bytes())
        expected_sequencing = _sequencing_rows(expected["sequencing"])
    except Exception as exc:
        record(
            "semantic:sequencing:read",
            False,
            str(exc),
            artifact="results/sequencing/sequencing_order.xlsx",
            field="read",
            expected="readable canonical sequencing workbook",
            actual=str(exc),
        )
    else:
        _compare_workbook_rows(
            record,
            artifact="results/sequencing/sequencing_order.xlsx",
            expected=expected_sequencing,
            actual=actual_sequencing,
        )
    summary_path = paths["report"] / "run_summary.json"
    try:
        actual_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except Exception as exc:
        record(
            "semantic:summary:read",
            False,
            str(exc),
            artifact="results/report/run_summary.json",
            field="read",
            expected="readable JSON summary",
            actual=str(exc),
        )
    else:
        record(
            "semantic:summary:content",
            actual_summary == expected["summary"],
            "run summary",
            artifact="results/report/run_summary.json",
            field="summary",
            expected=expected["summary"],
            actual=actual_summary,
        )


def audit_project(
    project_root: str | Path,
    *,
    migration_registry: str | Path = DEFAULT_MIGRATION_REGISTRY,
) -> dict[str, Any]:
    """Audit one project without creating, changing, or deleting any asset."""
    root = Path(project_root)
    checks: list[dict[str, str]] = []

    def record(
        check: str,
        passed: bool,
        detail: str,
        *,
        artifact: str = "",
        construct: str = "",
        field: str = "",
        expected: Any = "",
        actual: Any = "",
    ) -> None:
        reported_expected = expected if not passed else ""
        reported_actual = actual if not passed else ""
        checks.append(
            {
                "check": check,
                "status": "PASS" if passed else "FAIL",
                "detail": detail,
                "artifact": artifact,
                "construct": construct,
                "field": field,
                "expected": _report_value(reported_expected),
                "actual": _report_value(reported_actual),
            }
        )

    record("project_root", root.is_dir(), str(root))
    for relative in CANONICAL_DIRECTORIES:
        path = root / relative
        record(f"directory:{relative}", path.is_dir(), str(path))
    for relative in REQUIRED_PROJECT_FILES:
        path = root / relative
        record(f"file:{relative}", path.is_file(), str(path))

    manifest_path = canonical_project_paths(root)["manifest"]
    manifest = None
    if manifest_path.is_file():
        try:
            manifest = load_project_manifest(manifest_path)
        except ProjectManifestError as exc:
            record("manifest_schema", False, str(exc))
        else:
            project_id_matches = manifest["project_id"] == root.name
            record(
                "manifest_schema",
                True,
                f"contract_version={manifest['contract_version']}",
            )
            record(
                "manifest_project_id",
                project_id_matches,
                f"manifest={manifest['project_id']!r}; directory={root.name!r}",
            )
    else:
        record("manifest_schema", False, f"missing {manifest_path}")

    profile, registration, profile_error = _registered_validation_profile(
        root, migration_registry
    )
    record(
        "validation_profile",
        profile_error is None,
        profile if profile_error is None else profile_error,
        artifact="migration_registry",
        field="validation_profile",
        expected=profile,
        actual=profile if profile_error is None else profile_error,
    )
    if profile == LEGACY_COMPATIBILITY:
        record(
            "legacy_annotation_snapshot_provenance",
            bool(registration.get("annotation_snapshot_provenance")),
            str(registration.get("annotation_snapshot_provenance") or ""),
            artifact="migration_registry",
            field="annotation_snapshot_provenance",
            expected="registered immutable legacy annotation snapshot",
            actual=registration.get("annotation_snapshot_provenance"),
        )
    required_present = all((root / relative).is_file() for relative in REQUIRED_PROJECT_FILES)
    if manifest is not None and required_present and profile_error is None:
        _audit_semantics(root, manifest, profile, registration, record)

    failed = sum(item["status"] == "FAIL" for item in checks)
    return {
        "project": root.name,
        "project_root": str(root),
        "contract_version": manifest.get("contract_version") if manifest else None,
        "validation_profile": profile,
        "contract_generation": registration.get("contract_generation")
        if profile == LEGACY_COMPATIBILITY
        else "current",
        "status": "PASS" if failed == 0 else "FAIL",
        "failed_checks": failed,
        "checks": checks,
    }


def audit_primer_portfolio(
    migration_registry: str | Path = DEFAULT_MIGRATION_REGISTRY,
    *,
    skill_root: str | Path = DEFAULT_PRIMER_SKILL_ROOT,
) -> dict[str, Any]:
    """Validate every governed project; active migrations block completion."""
    checks: list[dict[str, str]] = []

    def record(check: str, passed: bool, detail: str) -> None:
        checks.append({"check": check, "status": "PASS" if passed else "FAIL", "detail": detail})

    migration = audit_migration_registry(migration_registry)
    record(
        "migration_registry",
        migration["status"] == "PASS",
        f"failed_checks={migration['failed_checks']}",
    )
    registry = load_migration_registry(migration_registry)
    root = Path(registry["project_root"])
    projects = []
    for entry in registry["projects"]:
        project_id = entry["project_id"]
        project_root = root / project_id
        manifest_path = project_root / "project.yaml"
        if manifest_path.is_file():
            audit = audit_project(project_root)
            status = audit["status"]
            detail = f"canonical contract; failed_checks={audit['failed_checks']}"
            mode = "CANONICAL"
        elif entry["migration_action"] in {"INDEX_ONLY", "KEEP"}:
            status = "PASS"
            detail = f"historical {entry['migration_action']}; source-of-truth indexed"
            mode = "HISTORICAL_INDEX"
        else:
            state = str(entry["migration_state"])
            status = "REVIEW" if "BLOCKED" in state or entry["migration_action"] == "REVIEW" else "PARTIAL"
            detail = f"{entry['migration_action']} required; state={entry['migration_state']}"
            mode = "MIGRATION_REQUIRED"
        record(f"project:{project_id}", status == "PASS", detail)
        projects.append(
            {
                "project_id": project_id,
                "project_status": entry["project_status"],
                "mode": mode,
                "status": status,
                "detail": detail,
            }
        )

    try:
        order_template = resolve_order_template()
    except Exception as exc:
        record("shared_order_template", False, str(exc))
    else:
        record("shared_order_template", True, order_template.template_name)
    try:
        sequencing_template = resolve_sequencing_template()
    except Exception as exc:
        record("shared_sequencing_template", False, str(exc))
    else:
        record(
            "shared_sequencing_template",
            True,
            f"{sequencing_template.source_path}; kind={sequencing_template.source_kind}; "
            f"sha256={sequencing_template.sha256}",
        )

    skill = Path(skill_root)
    required_skill_files = (skill / "SKILL.md", skill / "references/project-contract.md")
    missing_skill_files = [str(path) for path in required_skill_files if not path.is_file()]
    record(
        "primer_skill_contract",
        not missing_skill_files,
        "skill contract files present" if not missing_skill_files else str(missing_skill_files),
    )

    failed = sum(check["status"] == "FAIL" for check in checks)
    return {
        "status": "PASS" if failed == 0 else "FAIL",
        "normalization_complete": failed == 0,
        "failed_checks": failed,
        "projects": projects,
        "checks": checks,
        "migration_audit": migration,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit one Primer Project Contract directory")
    parser.add_argument("project_root", type=Path, nargs="?")
    parser.add_argument("--portfolio", action="store_true")
    parser.add_argument("--migration-registry", type=Path, default=DEFAULT_MIGRATION_REGISTRY)
    args = parser.parse_args()
    if args.portfolio:
        report = audit_primer_portfolio(args.migration_registry)
    elif args.project_root is not None:
        report = audit_project(args.project_root)
    else:
        parser.error("project_root is required unless --portfolio is used")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
