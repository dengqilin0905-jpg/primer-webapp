"""Read-only migration registry and audit for historical primer projects."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from project_contract import LEGACY_COMPATIBILITY, VALIDATION_PROFILES


from resource_config import resource_path
DEFAULT_MIGRATION_REGISTRY = resource_path('migration_registry',Path(__file__).resolve().parents[1]/'config/project_migrations.json')
ALLOWED_PROJECT_STATUSES = {"ACTIVE", "COMPLETED", "TEST_FIXTURE", "DUPLICATE_OBSOLETE"}
ALLOWED_MIGRATION_ACTIONS = {"KEEP", "MIGRATE", "INDEX_ONLY", "REVIEW"}
STATUS_ACTIONS = {
    "ACTIVE": {"MIGRATE"},
    "COMPLETED": {"KEEP", "INDEX_ONLY"},
    "TEST_FIXTURE": {"KEEP"},
    "DUPLICATE_OBSOLETE": {"REVIEW"},
}


class ProjectMigrationError(ValueError):
    """Raised when the migration registry is unsafe or incomplete."""


def load_migration_registry(path: str | Path = DEFAULT_MIGRATION_REGISTRY) -> dict[str, Any]:
    source = Path(path)
    try:
        data = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ProjectMigrationError(f"Unable to read migration registry {source}: {exc}") from exc
    if not isinstance(data, dict):
        raise ProjectMigrationError("Migration registry root must be an object")
    return data


def _safe_child(root: Path, relative: str) -> Path:
    candidate = (root / relative).resolve()
    try:
        candidate.relative_to(root.resolve())
    except ValueError as exc:
        raise ProjectMigrationError(f"Registered path escapes project root: {relative}") from exc
    return candidate


def audit_migration_registry(
    path: str | Path = DEFAULT_MIGRATION_REGISTRY,
) -> dict[str, Any]:
    """Audit registry coverage and source paths without modifying project files."""
    data = load_migration_registry(path)
    checks: list[dict[str, str]] = []

    def record(check: str, passed: bool, detail: str) -> None:
        checks.append({"check": check, "status": "PASS" if passed else "FAIL", "detail": detail})

    record("registry_version", data.get("registry_version") == 1, str(data.get("registry_version")))
    root = Path(str(data.get("project_root") or ""))
    record("project_root", root.is_dir(), str(root))
    projects = data.get("projects")
    if not isinstance(projects, list):
        projects = []
        record("projects_schema", False, "projects must be a list")
    else:
        record("projects_schema", True, f"projects={len(projects)}")

    excluded = data.get("excluded_entries")
    if not isinstance(excluded, list):
        excluded = []
        record("excluded_schema", False, "excluded_entries must be a list")
    excluded_names = {
        str(item.get("name")) for item in excluded if isinstance(item, dict) and item.get("name")
    }
    record("global_order_template_excluded", "order_primer.xlsx" in excluded_names, str(excluded_names))

    registered_ids = []
    project_reports = []
    for index, project in enumerate(projects):
        prefix = f"projects[{index}]"
        if not isinstance(project, dict):
            record(prefix, False, "project record must be an object")
            continue
        project_id = str(project.get("project_id") or "")
        status = str(project.get("project_status") or "")
        action = str(project.get("migration_action") or "")
        state = str(project.get("migration_state") or "")
        validation_profile = str(project.get("validation_profile") or "").strip()
        registered_ids.append(project_id)
        project_root = _safe_child(root, project_id) if project_id else root
        record(f"{project_id}:directory", bool(project_id) and project_root.is_dir(), str(project_root))
        record(f"{project_id}:status", status in ALLOWED_PROJECT_STATUSES, status)
        record(f"{project_id}:action", action in ALLOWED_MIGRATION_ACTIONS, action)
        record(
            f"{project_id}:classification",
            action in STATUS_ACTIONS.get(status, set()),
            f"{status}->{action}",
        )
        record(f"{project_id}:migration_state", bool(state), state)
        if validation_profile:
            record(
                f"{project_id}:validation_profile",
                validation_profile in VALIDATION_PROFILES,
                validation_profile,
            )
        if validation_profile == LEGACY_COMPATIBILITY:
            record(
                f"{project_id}:legacy_immutable",
                project.get("immutable") is True,
                str(project.get("immutable")),
            )
            record(
                f"{project_id}:contract_generation",
                bool(str(project.get("contract_generation") or "").strip()),
                str(project.get("contract_generation") or ""),
            )
            record(
                f"{project_id}:annotation_snapshot_provenance",
                bool(str(project.get("annotation_snapshot_provenance") or "").strip()),
                str(project.get("annotation_snapshot_provenance") or ""),
            )
        sources = project.get("source_of_truth")
        sources = sources if isinstance(sources, list) else []
        record(f"{project_id}:sources_declared", bool(sources), f"count={len(sources)}")
        missing_sources = []
        for relative in sources:
            source = _safe_child(project_root, str(relative))
            if not source.exists():
                missing_sources.append(str(relative))
        record(
            f"{project_id}:sources_exist",
            not missing_sources,
            "all source-of-truth paths exist" if not missing_sources else str(missing_sources),
        )
        issues = project.get("issues")
        issues_valid = isinstance(issues, list) and all(
            isinstance(issue, dict)
            and all(str(issue.get(field) or "").strip() for field in ("code", "detail", "disposition"))
            for issue in issues
        )
        record(f"{project_id}:issues", issues_valid, f"count={len(issues) if isinstance(issues, list) else 0}")
        project_reports.append(
            {
                "project_id": project_id,
                "project_status": status,
                "migration_action": action,
                "migration_state": state,
                "validation_profile": validation_profile or None,
                "immutable": project.get("immutable") is True,
                "contract_generation": project.get("contract_generation"),
                "contract_status": "MIGRATION_REQUIRED" if action == "MIGRATE" else "HISTORICAL_INDEX",
                "issue_codes": [issue.get("code") for issue in issues] if isinstance(issues, list) else [],
            }
        )

    duplicate_ids = sorted({item for item in registered_ids if registered_ids.count(item) > 1})
    record("unique_project_ids", not duplicate_ids and all(registered_ids), str(duplicate_ids))
    actual_directories = sorted(path.name for path in root.iterdir() if path.is_dir()) if root.is_dir() else []
    record(
        "all_project_directories_registered",
        sorted(registered_ids) == actual_directories,
        f"registered={sorted(registered_ids)}; actual={actual_directories}",
    )
    unexpected_entries = []
    if root.is_dir():
        for entry in root.iterdir():
            if entry.is_dir() or entry.name in excluded_names:
                continue
            unexpected_entries.append(entry.name)
    record("non_project_entries_excluded", not unexpected_entries, str(sorted(unexpected_entries)))

    failed = sum(check["status"] == "FAIL" for check in checks)
    return {
        "registry": str(Path(path)),
        "project_root": str(root),
        "status": "PASS" if failed == 0 else "FAIL",
        "failed_checks": failed,
        "projects": project_reports,
        "checks": checks,
    }
