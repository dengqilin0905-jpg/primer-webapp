"""Normalize shared presets and user backbones for the Design page."""

from __future__ import annotations

import io
import json
from pathlib import Path
from typing import Any

from backbone_presets import REPO_ROOT, get_backbone_preset
from seq_io import parse_uploaded_records


def _uploaded_bytes(path: Path):
    handle = io.BytesIO(path.read_bytes())
    handle.name = path.name
    return handle


def _load_shared_vector(preset: dict[str, Any]) -> dict[str, Any]:
    path = (REPO_ROOT / preset["vector"]["path"]).resolve()
    records = parse_uploaded_records(_uploaded_bytes(path))
    if len(records) != 1:
        raise ValueError(f"Expected one vector record for {preset['preset_id']}")
    return records[0]


def shared_preset_design_config(preset_id: str, *, version: int | None = None) -> dict[str, Any]:
    """Resolve and validate a registry preset into the engine-facing fields."""
    preset = get_backbone_preset(preset_id, version=version)
    vector = _load_shared_vector(preset)
    cloning = preset["cloning"]
    homology = cloning["homology"]
    fusion = cloning["fusion"]
    return {
        "backbone_source": "shared_preset",
        "backbone_name": preset["preset_id"],
        "backbone_label": f"{preset['preset_id']} [shared v{preset['preset_version']}]",
        "preset_id": preset["preset_id"],
        "preset_version": preset["preset_version"],
        "preset_status": preset["status"],
        "validation_level": preset["validation_level"],
        "approved_by_user": bool(preset.get("approved_by_user", False)),
        "vector_sha256": preset["vector"]["sha256"],
        "vector_sequence": vector["sequence"],
        "vector_features": vector.get("features") or [],
        "vector_features_json": json.dumps(vector.get("features") or [], ensure_ascii=False),
        "vector_topology": vector.get("topology") or preset["vector"].get("topology", ""),
        "vector_source_format": vector.get("source_format") or preset["vector"].get("source_format", ""),
        "vector_record_name": vector.get("record_name") or preset["vector"].get("record_name", ""),
        "insertion_mode": cloning["insertion_mode"],
        "restriction_site_f": cloning["restriction_site_f"]["sequence"],
        "restriction_site_r": cloning["restriction_site_r"]["sequence"],
        "enzyme_f": cloning["restriction_site_f"]["enzyme"],
        "enzyme_r": cloning["restriction_site_r"]["enzyme"],
        "enzyme_sites": f'{cloning["restriction_site_f"]["enzyme"]}/{cloning["restriction_site_r"]["enzyme"]}',
        "homo_f": homology["homo_f"],
        "homo_r": homology["homo_r"],
        "homology_arm_length": homology["arm_length"],
        "junction_extra_bases_f": 0,
        "junction_extra_bases_r": 0,
        "vector_map_note": preset["provenance"].get("validation_note", ""),
        "construct_type": fusion["construct_type"],
        "fusion_order": fusion["order"],
        "fusion_architecture": fusion["architecture"],
        "remove_stop_default": bool(preset["design_defaults"].get("remove_stop", False)),
        "primer_engine_frame_check": bool(preset["qc"].get("primer_engine_frame_check", False)),
        "fusion_frame_validation": bool(preset["qc"].get("fusion_frame_validation", False)),
        "validated_replacement": dict(cloning.get("replacement") or {}),
        "provenance": {
            "source_type": "shared_preset",
            "preset_id": preset["preset_id"],
            "preset_version": preset["preset_version"],
            "preset_status": preset["status"],
            "validation_level": preset["validation_level"],
            "approved_by_user": bool(preset.get("approved_by_user", False)),
            "vector_sha256": preset["vector"]["sha256"],
        },
    }


def resolve_batch_backbone(
    preset_id: str, *, version: int | None = None, overrides: dict[str, Any] | None = None
) -> dict[str, Any]:
    """Resolve a shared preset for future Batch/Codex execution.

    The vector identity and validated cloning profile are immutable here.  A
    construct or project may explicitly override only the CDS stop rule; a
    changed vector, site, arm, boundary, or fusion architecture requires a
    different preset version or a custom backbone.
    """
    requested = dict(overrides or {})
    allowed = {"remove_stop_override"}
    forbidden = sorted(set(requested) - allowed)
    if forbidden:
        raise ValueError(
            "Shared preset cloning fields cannot be overridden in Batch/Codex; "
            f"use another preset version or a custom backbone (got: {', '.join(forbidden)})"
        )
    override = requested.get("remove_stop_override")
    if override is not None and not isinstance(override, bool):
        raise ValueError("remove_stop_override must be True, False, or omitted")
    config = dict(shared_preset_design_config(preset_id, version=version))
    effective_stop, differs = effective_remove_stop(config, "CDS", override)
    config["effective_remove_stop"] = effective_stop
    config["remove_stop_overridden"] = differs
    return config


def custom_backbone_design_config(row) -> dict[str, Any]:
    """Resolve one selected user backbone; never use another row's defaults."""
    data = dict(row)
    try:
        defaults = json.loads(data.get("defaults_json") or "{}")
    except (TypeError, json.JSONDecodeError):
        defaults = {}
    return {
        "backbone_source": "user_backbone",
        "backbone_id": int(data["id"]),
        "backbone_name": data["name"],
        "backbone_label": data.get("display") or data["name"],
        "preset_id": None,
        "preset_version": None,
        "preset_status": None,
        "validation_level": None,
        "approved_by_user": False,
        "vector_sha256": None,
        "vector_sequence": data.get("vector_sequence") or "",
        "vector_features": json.loads(data.get("vector_features_json") or "[]"),
        "vector_features_json": data.get("vector_features_json") or "[]",
        "vector_topology": data.get("vector_topology") or "",
        "vector_source_format": data.get("vector_source_format") or "",
        "vector_record_name": data.get("vector_record_name") or "",
        "insertion_mode": data.get("insertion_mode") or "linear_join",
        "restriction_site_f": data.get("restriction_site_f") or "",
        "restriction_site_r": data.get("restriction_site_r") or "",
        "enzyme_f": (data.get("enzyme_sites") or "").split("/", 1)[0],
        "enzyme_r": (data.get("enzyme_sites") or "").split("/", 1)[-1],
        "enzyme_sites": data.get("enzyme_sites") or "",
        "homo_f": data.get("homo_f") or "",
        "homo_r": data.get("homo_r") or "",
        "homology_arm_length": len(data.get("homo_f") or "") or None,
        "junction_extra_bases_f": int(data.get("junction_extra_bases_f") or 0),
        "junction_extra_bases_r": int(data.get("junction_extra_bases_r") or 0),
        "vector_map_note": data.get("vector_map_note") or "",
        "construct_type": defaults.get("construct_type", "CDS_no_tag"),
        "fusion_order": defaults.get("fusion_order"),
        "fusion_architecture": defaults.get("fusion_architecture"),
        "remove_stop_default": bool(defaults.get("remove_stop", True)),
        "primer_engine_frame_check": bool(defaults.get("enable_frame_check", False)),
        "fusion_frame_validation": bool(defaults.get("fusion_frame_validation", False)),
        "validated_replacement": {},
        "design_defaults": defaults,
        "provenance": {"source_type": "user_backbone", "backbone_id": int(data["id"])},
    }


def effective_remove_stop(config: dict[str, Any], seq_type: str, override: bool | None = None) -> tuple[bool, bool]:
    """Return (effective value, whether it differs from the selected default)."""
    if str(seq_type).upper() != "CDS":
        return False, False
    default = bool(config.get("remove_stop_default", False))
    value = default if override is None else bool(override)
    return value, value != default


def config_provenance(config: dict[str, Any], effective_stop: bool) -> dict[str, Any]:
    return {
        "backbone_source": config["backbone_source"],
        "preset_id": config.get("preset_id"),
        "preset_version": config.get("preset_version"),
        "preset_status": config.get("preset_status"),
        "validation_level": config.get("validation_level"),
        "approved_by_user": config.get("approved_by_user"),
        "vector_sha256": config.get("vector_sha256"),
        "effective_remove_stop": bool(effective_stop),
        "validated_homology_arm_length": config.get("homology_arm_length"),
        "fusion_order": config.get("fusion_order"),
        "fusion_architecture": config.get("fusion_architecture"),
    }
