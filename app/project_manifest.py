"""Read and validate a Primer Project Contract manifest.

Checked-in and generated ``project.yaml`` files use the JSON-compatible subset
of YAML, matching the existing backbone registry convention and avoiding a new
runtime dependency.
"""

from __future__ import annotations

import copy
import json
import os
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from project_contract import (
    ALLOWED_OVERRIDE_TYPES,
    CANONICAL_OUTPUTS,
    CONTRACT_VERSION,
    canonical_species_abbreviation,
)
from naming import build_operational_label


class ProjectManifestError(ValueError):
    """Raised when a project manifest cannot satisfy the contract."""


class ProjectManifestConflictError(ProjectManifestError):
    """Raised when a supplied manifest conflicts with an existing control file."""


@dataclass(frozen=True)
class ProjectManifestResolution:
    manifest: dict[str, Any]
    path: Path
    action: str


REQUIRED_CONSTRUCT_FIELDS = (
    "gene_id",
    "symbol",
    "canonical_identity",
    "operational_label",
    "full_construct_name",
    "backbone",
    "preset_id",
    "preset_version",
    "remove_stop",
    "physical_plasmid_status",
    "overrides",
    "status",
)


def _nonempty(value: Any) -> bool:
    return value not in (None, "", [], {})


def _validate_created_at(value: Any, errors: list[str]) -> None:
    if not isinstance(value, str) or not value.strip():
        errors.append("created_at must be a non-empty ISO-8601 string")
        return
    try:
        datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        errors.append("created_at must be ISO-8601 compatible")


def manifest_errors(data: Any) -> list[str]:
    """Return deterministic validation errors without raising or writing."""
    if not isinstance(data, dict):
        return ["manifest root must be an object"]

    errors: list[str] = []
    if data.get("contract_version") != CONTRACT_VERSION:
        errors.append(f"contract_version must be {CONTRACT_VERSION}")
    for field in ("project_id", "project_type"):
        if not _nonempty(data.get(field)):
            errors.append(f"missing required field: {field}")
    _validate_created_at(data.get("created_at"), errors)

    species = data.get("species")
    species_values = species if isinstance(species, list) else [species]
    if not species_values or any(not _nonempty(item) for item in species_values):
        errors.append("species must contain at least one value")
    else:
        for item in species_values:
            try:
                canonical_species_abbreviation(str(item))
            except ValueError as exc:
                errors.append(str(exc))

    constructs = data.get("constructs")
    if not isinstance(constructs, list):
        errors.append("constructs must be a list")
    else:
        for index, construct in enumerate(constructs):
            prefix = f"constructs[{index}]"
            if not isinstance(construct, dict):
                errors.append(f"{prefix} must be an object")
                continue
            for field in REQUIRED_CONSTRUCT_FIELDS:
                if field not in construct:
                    errors.append(f"{prefix} missing required field: {field}")
            if "remove_stop" in construct and not isinstance(construct["remove_stop"], bool):
                errors.append(f"{prefix}.remove_stop must be boolean")
            overrides = construct.get("overrides")
            if not isinstance(overrides, list):
                errors.append(f"{prefix}.overrides must be a list")
            else:
                for override_index, override in enumerate(overrides):
                    oprefix = f"{prefix}.overrides[{override_index}]"
                    if not isinstance(override, dict):
                        errors.append(f"{oprefix} must be an object")
                        continue
                    for field in ("override_type", "override_reason", "override_value"):
                        if not _nonempty(override.get(field)):
                            errors.append(f"{oprefix} missing required field: {field}")
                    override_type = override.get("override_type")
                    if override_type and override_type not in ALLOWED_OVERRIDE_TYPES:
                        errors.append(f"{oprefix}.override_type is not allowed: {override_type}")

            construct_species = construct.get("species")
            if not construct_species and len(species_values) == 1:
                construct_species = species_values[0]
            if construct_species and all(field in construct for field in ("symbol", "backbone")):
                try:
                    expected_label = build_operational_label(
                        construct.get("symbol"),
                        species=str(construct_species),
                        identifier=construct.get("canonical_identity") or None,
                        backbone=construct.get("backbone"),
                    )
                except ValueError as exc:
                    errors.append(f"{prefix}.operational_label cannot be resolved: {exc}")
                else:
                    if construct.get("operational_label") != expected_label:
                        errors.append(f"{prefix}.operational_label must be {expected_label!r}")

    outputs = data.get("outputs")
    if not isinstance(outputs, dict):
        errors.append("outputs must be an object")
    else:
        for key, expected in CANONICAL_OUTPUTS.items():
            if outputs.get(key) != expected:
                errors.append(f"outputs.{key} must be {expected!r}")

    provenance = data.get("provenance")
    if not isinstance(provenance, dict):
        errors.append("provenance must be an object")
    else:
        for field in ("identity_sources", "vector_hash", "template_source"):
            if field not in provenance:
                errors.append(f"provenance missing required field: {field}")
        if "identity_sources" in provenance and not isinstance(provenance["identity_sources"], list):
            errors.append("provenance.identity_sources must be a list")

    return errors


def validate_project_manifest(data: Any) -> dict[str, Any]:
    """Validate and return a defensive manifest copy."""
    errors = manifest_errors(data)
    if errors:
        raise ProjectManifestError("; ".join(errors))
    return copy.deepcopy(data)


def _strip_yaml_comments(text: str) -> str:
    """Strip comments while preserving hash characters inside JSON strings."""
    output = []
    for line in text.splitlines():
        in_string = False
        escaped = False
        kept = []
        for character in line:
            if escaped:
                kept.append(character)
                escaped = False
                continue
            if character == "\\" and in_string:
                kept.append(character)
                escaped = True
                continue
            if character == '"':
                kept.append(character)
                in_string = not in_string
                continue
            if character == "#" and not in_string:
                break
            kept.append(character)
        output.append("".join(kept))
    return "\n".join(output)


def load_project_manifest_bytes(data: bytes, *, source: str = "project.yaml") -> dict[str, Any]:
    """Load the documented JSON-compatible YAML subset, including comments."""
    try:
        value = json.loads(_strip_yaml_comments(data.decode("utf-8")))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ProjectManifestError(f"Unable to read project manifest {source}: {exc}") from exc
    return validate_project_manifest(value)


def load_project_manifest(path: str | Path) -> dict[str, Any]:
    """Load a JSON-compatible YAML manifest and validate it."""
    source = Path(path)
    try:
        data = source.read_bytes()
    except OSError as exc:
        raise ProjectManifestError(f"Unable to read project manifest {source}: {exc}") from exc
    return load_project_manifest_bytes(data, source=str(source))


def _canonical_manifest_bytes(manifest: dict[str, Any]) -> bytes:
    return json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True).encode("utf-8") + b"\n"


def _atomic_create_manifest(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=".project.yaml.", dir=str(path.parent))
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError as exc:
            raise ProjectManifestConflictError(
                f"Project manifest appeared during initialization: {path}"
            ) from exc
    finally:
        temporary.unlink(missing_ok=True)


def resolve_project_manifest(
    project_root: str | Path,
    supplied_manifest: dict[str, Any] | None = None,
) -> ProjectManifestResolution:
    """Create, reuse, or reject the control manifest without rewriting it."""
    root = Path(project_root)
    path = root / "project.yaml"
    supplied = validate_project_manifest(supplied_manifest) if supplied_manifest is not None else None
    if path.exists():
        existing = load_project_manifest(path)
        if supplied is not None and existing != supplied:
            raise ProjectManifestConflictError(
                f"Supplied manifest conflicts with existing project manifest: {path}"
            )
        return ProjectManifestResolution(existing, path, "reused")
    if supplied is None:
        raise ProjectManifestError(f"Project manifest does not exist and none was supplied: {path}")
    _atomic_create_manifest(path, _canonical_manifest_bytes(supplied))
    return ProjectManifestResolution(supplied, path, "created")
