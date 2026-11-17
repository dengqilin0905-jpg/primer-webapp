"""Read-only TF/RTF identity and physical-plasmid preflight.

Identity assignment and physical build state are intentionally independent:
designs, maps, database rows, and ChIP data can prove an identity, but they do
not prove that a usable physical plasmid exists.
"""

from __future__ import annotations

import csv
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from resource_config import resource_path,REPO_ROOT
from typing import Iterable, Mapping

from project_contract import (
    HARD_REQUIRED_IDENTITY_SOURCE_KINDS,
    IDENTITY_PREFLIGHT_SOURCE_KINDS,
    IDENTITY_SOURCE_SEVERITY,
    PROJECTS_ROOT,
)


RICE_IDENTITY_REGISTRY = resource_path('identity_registry',REPO_ROOT/'resources/examples/identity.tsv')
CANONICAL_CONSTRUCT_MASTER = resource_path('construct_master',REPO_ROOT/'resources/examples/identity.tsv')

IDENTITY_COLUMNS = (
    "canonical_identity",
    "current_id",
    "current_ID",
    "construct_ID",
    "identity",
    "identifier",
)
GENE_ID_COLUMNS = ("gene_ID", "gene_id", "Gene_ID")
PHYSICAL_STATUS_COLUMNS = (
    "physical_plasmid_status",
    "existing_physical_tagging_plasmid",
    "physical_status",
)

AVAILABLE_PHYSICAL_STATES = {"available", "confirmed", "yes", "pass", "built", "verified", "exists"}
MISSING_PHYSICAL_STATES = {"missing", "not_found", "failed", "no", "unavailable", "lost"}
UNCERTAIN_PHYSICAL_STATES = {"", "na", "n/a", "uncertain", "unknown", "not_checked", "pending"}

_TF_ID_RE = re.compile(r"^(?:Os)?TF(\d+)$", re.IGNORECASE)
_IDENTITY_KEY_RE = re.compile(r"^(?:Os)?(RTF|TF)(\d+)$", re.IGNORECASE)


class IdentityAuditError(ValueError):
    """Raised when identity evidence cannot be read safely."""


@dataclass(frozen=True)
class IdentityEvidence:
    gene_id: str
    canonical_identity: str
    source_kind: str
    source_path: str
    physical_plasmid_status: str = "unknown"


def _first(row: Mapping[str, str], columns: Iterable[str]) -> str:
    for column in columns:
        value = str(row.get(column) or "").strip()
        if value and value.upper() != "NA":
            return value
    return ""


def identity_key(identity: str) -> str:
    """Normalize species-prefixed and historical TF/RTF spellings for dedupe."""
    value = str(identity or "").strip()
    match = _IDENTITY_KEY_RE.fullmatch(value)
    if not match:
        return value.upper()
    return f"{match.group(1).upper()}{int(match.group(2))}"


def resolve_identity_column(columns: Iterable[str]) -> tuple[str, ...]:
    """Return exact identity aliases in canonical priority order.

    Returning every present alias is intentional: row loading must compare
    their contents and may not silently select one conflicting field.
    """
    available = set(columns)
    resolved = tuple(column for column in IDENTITY_COLUMNS if column in available)
    if not resolved:
        raise IdentityAuditError(
            f"Identity source lacks an identity column; accepted aliases: {', '.join(IDENTITY_COLUMNS)}"
        )
    return resolved


def _identity_from_row(
    row: Mapping[str, str],
    identity_columns: Iterable[str],
    *,
    source: Path,
    row_number: int,
) -> str:
    values = [
        (column, str(row.get(column) or "").strip())
        for column in identity_columns
        if str(row.get(column) or "").strip()
        and str(row.get(column) or "").strip().upper() != "NA"
    ]
    keys = {identity_key(value) for _, value in values}
    if len(keys) > 1:
        detail = ", ".join(f"{column}={value!r}" for column, value in values)
        raise IdentityAuditError(f"Conflicting identity aliases in {source} row {row_number}: {detail}")
    return values[0][1] if values else ""


def classify_physical_status(value: str | None) -> str:
    """Collapse free-text inventory states without treating design as physical."""
    normalized = str(value or "").strip().lower().replace(" ", "_").replace("-", "_")
    if normalized in AVAILABLE_PHYSICAL_STATES:
        return "available"
    if normalized in MISSING_PHYSICAL_STATES:
        return "missing"
    if normalized in UNCERTAIN_PHYSICAL_STATES:
        return "uncertain"
    return "uncertain"


def load_identity_tsv(path: str | Path, *, source_kind: str) -> list[IdentityEvidence]:
    """Load identity evidence from a TSV using known canonical column aliases."""
    source = Path(path)
    if source_kind not in IDENTITY_PREFLIGHT_SOURCE_KINDS:
        raise IdentityAuditError(f"Unsupported identity source kind: {source_kind}")
    try:
        handle = source.open("r", encoding="utf-8-sig", newline="")
    except OSError as exc:
        raise IdentityAuditError(f"Unable to read identity source {source}: {exc}") from exc
    evidence: list[IdentityEvidence] = []
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise IdentityAuditError(f"Identity source has no TSV header: {source}")
        if not any(column in reader.fieldnames for column in GENE_ID_COLUMNS):
            raise IdentityAuditError(f"Identity source lacks a gene ID column: {source}")
        try:
            identity_columns = resolve_identity_column(reader.fieldnames)
        except IdentityAuditError as exc:
            raise IdentityAuditError(f"{exc}: {source}") from exc
        for row_number, row in enumerate(reader, start=2):
            gene_id = _first(row, GENE_ID_COLUMNS)
            identity = _identity_from_row(
                row,
                identity_columns,
                source=source,
                row_number=row_number,
            )
            if not gene_id or not identity:
                continue
            evidence.append(
                IdentityEvidence(
                    gene_id=gene_id,
                    canonical_identity=identity,
                    source_kind=source_kind,
                    source_path=str(source),
                    physical_plasmid_status=_first(row, PHYSICAL_STATUS_COLUMNS) or "unknown",
                )
            )
    return evidence


def collect_identity_evidence(
    sources: Mapping[str, Iterable[str | Path]],
) -> list[IdentityEvidence]:
    """Load explicitly declared sources; the function never searches implicitly."""
    evidence: list[IdentityEvidence] = []
    for source_kind, paths in sources.items():
        for path in paths:
            evidence.extend(load_identity_tsv(path, source_kind=source_kind))
    return evidence


def discover_active_project_evidence(
    *,
    projects_root: str | Path = PROJECTS_ROOT,
    exclude_project_root: str | Path | None = None,
) -> list[IdentityEvidence]:
    """Read active-project identities without requiring a user-supplied plan.

    Registered ACTIVE legacy projects contribute their declared construct plans.
    Canonical projects contribute their manifest constructs.  Missing or invalid
    supplementary plans are ignored here; their absence is reported by source
    status and cannot downgrade the hard 2026 identity registry.
    """
    root = Path(projects_root)
    excluded = Path(exclude_project_root).resolve() if exclude_project_root else None
    evidence: list[IdentityEvidence] = []
    seen_paths: set[Path] = set()
    registry = root.parent / "repo" / "config" / "project_migrations.json"
    try:
        records = json.loads(registry.read_text(encoding="utf-8")).get("projects", [])
    except (OSError, json.JSONDecodeError, AttributeError):
        records = []
    for record in records:
        if not isinstance(record, dict) or record.get("project_status") != "ACTIVE":
            continue
        project = root / str(record.get("project_id") or "")
        if excluded and project.resolve() == excluded:
            continue
        for relative in record.get("source_of_truth", []):
            path = project / str(relative)
            if path.name != "construct_plan.tsv" or not path.is_file() or path in seen_paths:
                continue
            try:
                evidence.extend(load_identity_tsv(path, source_kind="active_construct_plan"))
            except IdentityAuditError:
                continue
            seen_paths.add(path)
    if not root.is_dir():
        return evidence
    for manifest in root.glob("*/project.yaml"):
        if excluded and manifest.parent.resolve() == excluded:
            continue
        try:
            data = json.loads(manifest.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        for construct in data.get("constructs", []) if isinstance(data, dict) else []:
            if not isinstance(construct, dict):
                continue
            gene_id = str(construct.get("gene_id") or "").strip()
            identity = str(construct.get("canonical_identity") or "").strip()
            if gene_id and identity:
                evidence.append(IdentityEvidence(
                    gene_id=gene_id,
                    canonical_identity=identity,
                    source_kind="active_construct_plan",
                    source_path=str(manifest),
                    physical_plasmid_status=str(construct.get("physical_plasmid_status") or "unknown"),
                ))
    return evidence


def next_available_tf_id(evidence: Iterable[IdentityEvidence]) -> str:
    """Calculate, but never persist, the next TF number across all evidence."""
    numbers = []
    for item in evidence:
        match = _TF_ID_RE.fullmatch(item.canonical_identity.strip())
        if match:
            numbers.append(int(match.group(1)))
    return f"TF{max(numbers, default=0) + 1}"


def reserved_identities(evidence: Iterable[IdentityEvidence]) -> set[str]:
    """Return canonical spellings reserved by any evidence, even on BLOCK."""
    priority = {kind: index for index, kind in enumerate(IDENTITY_PREFLIGHT_SOURCE_KINDS)}
    reserved: dict[str, IdentityEvidence] = {}
    for item in evidence:
        key = identity_key(item.canonical_identity)
        current = reserved.get(key)
        if current is None or priority.get(item.source_kind, 999) < priority.get(current.source_kind, 999):
            reserved[key] = item
    return {item.canonical_identity for item in reserved.values()}


def _preferred_identity(matches: list[IdentityEvidence]) -> str:
    priority = {kind: index for index, kind in enumerate(IDENTITY_PREFLIGHT_SOURCE_KINDS)}
    return min(matches, key=lambda item: priority.get(item.source_kind, 999)).canonical_identity


def audit_identity(
    gene_id: str,
    evidence: Iterable[IdentityEvidence],
    *,
    requested_identity: str | None = None,
    required_source_kinds: Iterable[str] = HARD_REQUIRED_IDENTITY_SOURCE_KINDS,
    respect_existing_physical: bool = False,
) -> dict:
    """Return a fail-closed decision before allocation, design, or export."""
    records = list(evidence)
    present_source_kinds = {item.source_kind for item in records}
    missing_source_kinds = [kind for kind in required_source_kinds if kind not in present_source_kinds]
    matches = [item for item in records if item.gene_id == str(gene_id).strip()]
    identity_keys = {identity_key(item.canonical_identity) for item in matches}
    requested_key = identity_key(requested_identity or "")
    existing_identity = _preferred_identity(matches) if matches else ""

    decision = "READY_FOR_NEW_IDENTITY"
    reason = "No canonical identity was found in the declared evidence sources"
    physical_status = "uncertain"

    if missing_source_kinds:
        decision = "BLOCK_INCOMPLETE_SOURCES"
        reason = "Hard-required identity preflight sources are missing"
    elif len(identity_keys) > 1:
        decision = "BLOCK_IDENTITY_CONFLICT"
        reason = f"Gene maps to conflicting canonical identities: {sorted(identity_keys)}"
    elif existing_identity and requested_identity and requested_key not in identity_keys:
        decision = "BLOCK_IDENTITY_CONFLICT"
        reason = f"Requested identity {requested_identity!r} conflicts with existing {existing_identity!r}"
    elif existing_identity:
        inventory_states = {
            classify_physical_status(item.physical_plasmid_status)
            for item in matches
            if item.source_kind == "physical_inventory"
        }
        considered_states = inventory_states or {
            classify_physical_status(item.physical_plasmid_status) for item in matches
        }
        decisive_states = considered_states - {"uncertain"}
        if len(decisive_states) > 1 and respect_existing_physical:
            decision = "BLOCK_PHYSICAL_CONFLICT"
            reason = "Physical inventory contains conflicting available/missing states"
        elif "available" in decisive_states and respect_existing_physical:
            physical_status = "available"
            decision = "SKIP_REBUILD"
            reason = "Canonical identity and usable physical plasmid already exist"
        elif "missing" in decisive_states:
            physical_status = "missing"
            decision = "REBUILD_EXISTING_IDENTITY"
            reason = "Physical plasmid is missing/failed; rebuild must retain the existing identity"
        else:
            physical_status = "available" if "available" in decisive_states else "unknown"
            decision = "REBUILD_EXISTING_IDENTITY"
            reason = "Canonical identity retained; physical plasmid status does not block an explicit rebuild"

    present = {item.source_kind for item in records}
    source_status = {
        kind: (
            "AVAILABLE" if kind in present
            else "MISSING_HARD_REQUIRED" if IDENTITY_SOURCE_SEVERITY[kind] == "HARD_REQUIRED"
            else "NOT_FOUND" if IDENTITY_SOURCE_SEVERITY[kind] == "SUPPLEMENTARY"
            else "NOT_VERIFIED"
        )
        for kind in IDENTITY_PREFLIGHT_SOURCE_KINDS
    }

    return {
        "gene_id": str(gene_id).strip(),
        "requested_identity": requested_identity or "",
        "canonical_identity": existing_identity,
        "identity_key": identity_key(existing_identity),
        "physical_plasmid_status": physical_status,
        "decision": decision,
        "reason": reason,
        "missing_source_kinds": missing_source_kinds,
        "source_status": source_status,
        "next_available_tf_id": next_available_tf_id(records),
        "reserved_identities": sorted(reserved_identities(records)),
        "evidence": [asdict(item) for item in matches],
    }
