"""Pure-Python loader and validator for the shared backbone preset registry."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Any

from seq_io import parse_uploaded_records


from resource_config import resource_path
SOFTWARE_ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = resource_path('backbone_registry',SOFTWARE_ROOT/'config/backbone_presets.yaml')
REPO_ROOT = resource_path('backbone_root',REGISTRY_PATH.parent.parent)


class BackbonePresetError(ValueError):
    """Raised when a preset is missing or fails source validation."""


def _read_registry(path: Path = REGISTRY_PATH) -> dict[str, Any]:
    try:
        raw = path.read_text(encoding="utf-8")
        # The checked-in .yaml is JSON-compatible YAML so production loading
        # does not require PyYAML, which is not an application dependency.
        data = json.loads(raw)
    except (OSError, json.JSONDecodeError) as exc:
        raise BackbonePresetError(f"Unable to read backbone registry: {path}: {exc}") from exc
    if not isinstance(data, dict) or not isinstance(data.get("presets"), list):
        raise BackbonePresetError("Backbone registry must contain a presets list")
    return data


def load_backbone_registry(*, validate: bool = True, path: Path = REGISTRY_PATH) -> dict[str, Any]:
    """Load the complete registry, optionally validating every preset."""
    data = _read_registry(path)
    if validate:
        for preset in data["presets"]:
            validate_backbone_preset(preset, registry_path=path)
    return copy.deepcopy(data)


def list_backbone_presets(*, path: Path = REGISTRY_PATH) -> list[dict[str, Any]]:
    """Return one current validated version for each production preset."""
    data = load_backbone_registry(path=path)
    current_versions = data.get("current_versions")
    if not isinstance(current_versions, dict):
        raise BackbonePresetError("Backbone registry must define current_versions")
    return [
        get_backbone_preset(preset_id, version=version, path=path)
        for preset_id, version in current_versions.items()
    ]


def get_backbone_preset(
    preset_id: str, *, version: int | None = None, path: Path = REGISTRY_PATH
) -> dict[str, Any]:
    data = load_backbone_registry(path=path)
    query = str(preset_id).strip()
    matches = [p for p in data["presets"] if p.get("preset_id") == query or query in p.get("aliases", [])]
    if version is not None:
        matches = [p for p in matches if p.get("preset_version") == version]
    else:
        preset_ids = {p.get("preset_id") for p in matches}
        if len(preset_ids) != 1:
            raise BackbonePresetError(f"Ambiguous backbone preset: {preset_id!r}")
        canonical_id = preset_ids.pop()
        current_versions = data.get("current_versions")
        if not isinstance(current_versions, dict) or canonical_id not in current_versions:
            raise BackbonePresetError(f"No current version configured for backbone preset: {canonical_id!r}")
        current_version = current_versions[canonical_id]
        matches = [p for p in matches if p.get("preset_version") == current_version]
    if not matches:
        raise BackbonePresetError(f"Backbone preset not found: {preset_id!r} version={version!r}")
    if len(matches) > 1:
        raise BackbonePresetError(f"Ambiguous backbone preset: {preset_id!r} version={version!r}")
    return matches[0]


def _required(preset: dict[str, Any], key: str) -> Any:
    value = preset.get(key)
    if value in (None, "", []):
        raise BackbonePresetError(f"Missing required preset field: {key}")
    return value


def _sequence_from_asset(path: Path) -> tuple[str, dict[str, Any]]:
    class UploadedBytes:
        def __init__(self, content: bytes, name: str):
            self._content, self.name = content, name

        def read(self) -> bytes:
            return self._content

    records = parse_uploaded_records(UploadedBytes(path.read_bytes(), path.name))
    if len(records) != 1:
        raise BackbonePresetError(f"Expected one vector record in {path}, found {len(records)}")
    return records[0]["sequence"].upper(), records[0]


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def validate_backbone_preset(preset: dict[str, Any], *, registry_path: Path = REGISTRY_PATH) -> bool:
    """Validate identity, vector content, sites, boundary, and cloning profile."""
    for key in ("preset_id", "preset_version", "status", "vector", "cloning", "design_defaults", "qc"):
        _required(preset, key)
    status = preset["status"]
    if status not in {"historical_validated", "current", "validated"}:
        raise BackbonePresetError(f"Unsupported preset status: {preset.get('preset_id')}={status!r}")
    validation_level = _required(preset, "validation_level")
    if status == "historical_validated" and validation_level != "historical":
        raise BackbonePresetError(f"Historical preset must declare validation_level=historical: {preset.get('preset_id')}")
    if status == "current":
        if validation_level != "computational" or preset.get("approved_by_user") is not True:
            raise BackbonePresetError(
                f"Current preset must be user-approved computational profile: {preset.get('preset_id')}"
            )
    if status == "validated" and validation_level != "validated":
        raise BackbonePresetError(f"Validated preset has invalid validation level: {preset.get('preset_id')}")
    vector = preset["vector"]
    cloning = preset["cloning"]
    for key in ("path", "sha256", "sequence_length", "topology", "source_format"):
        _required(vector, key)
    asset_root = REPO_ROOT if Path(registry_path).resolve()==REGISTRY_PATH.resolve() else Path(registry_path).parent.parent
    asset = (asset_root / vector["path"]).resolve()
    if not asset.is_file():
        raise BackbonePresetError(f"Vector asset does not exist: {asset}")
    sequence, record = _sequence_from_asset(asset)
    digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    if digest != vector["sha256"]:
        raise BackbonePresetError(f"Vector sequence hash mismatch for {preset['preset_id']}: {digest}")
    if len(sequence) != vector["sequence_length"]:
        raise BackbonePresetError(f"Vector length mismatch for {preset['preset_id']}")
    if record.get("topology") and record["topology"] != vector["topology"]:
        raise BackbonePresetError(f"Vector topology mismatch for {preset['preset_id']}")
    homology = cloning.get("homology", {})
    arm_length = homology.get("arm_length")
    homo_f, homo_r = homology.get("homo_f", ""), homology.get("homo_r", "")
    if not isinstance(arm_length, int) or len(homo_f) != arm_length or len(homo_r) != arm_length:
        raise BackbonePresetError(f"Invalid homology profile for {preset['preset_id']}")
    replacement = cloning.get("replacement", {})
    start, end = replacement.get("boundary_start"), replacement.get("boundary_end")
    if not isinstance(start, int) or not isinstance(end, int) or not (0 <= start < end <= len(sequence)):
        raise BackbonePresetError(f"Invalid replacement boundary for {preset['preset_id']}")
    left = cloning.get("restriction_site_f", {}).get("sequence", "").upper()
    right = cloning.get("restriction_site_r", {}).get("sequence", "").upper()
    if not left or not right or sequence[start - len(left):start] != left or sequence[end:end + len(right)] != right:
        raise BackbonePresetError(f"Restriction sites do not flank replacement boundary for {preset['preset_id']}")
    if sequence[start - arm_length:start] != homo_f:
        raise BackbonePresetError(f"Forward homology arm mismatch for {preset['preset_id']}")
    if _revcomp(sequence[end:end + arm_length]) != homo_r:
        raise BackbonePresetError(f"Reverse homology arm mismatch for {preset['preset_id']}")
    if not cloning.get("fusion", {}).get("order") or not cloning.get("fusion", {}).get("construct_type"):
        raise BackbonePresetError(f"Missing fusion profile for {preset['preset_id']}")
    return True
