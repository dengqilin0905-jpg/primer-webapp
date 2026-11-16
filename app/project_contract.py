"""Canonical Primer Project Contract constants and path helpers.

This module is deliberately side-effect free.  It defines the filesystem and
table contracts shared by Web and batch orchestration without creating or
modifying a project.
"""

from __future__ import annotations

from pathlib import Path
from resource_config import resource_path,state_root,REPO_ROOT


CONTRACT_VERSION = 1
CANONICAL_RESULT_VERSION = 1
CANONICAL_NORMALIZATION_VERSION = 1
CURRENT_STRICT = "CURRENT_STRICT"
LEGACY_COMPATIBILITY = "LEGACY_COMPATIBILITY"
VALIDATION_PROFILES = (CURRENT_STRICT, LEGACY_COMPATIBILITY)
CANONICAL_RESULT_SECTIONS = (
    "identity",
    "annotation",
    "naming",
    "sequence_and_backbone_provenance",
    "primers",
    "lengths",
    "rescue",
    "map_input",
    "order_input",
    "sequencing_input",
    "warnings",
    "provenance",
    "qc_rows",
)
PROJECTS_ROOT = resource_path('projects_root',state_root()/'projects')
GLOBAL_ORDER_TEMPLATE = resource_path('order_template',REPO_ROOT/'resources/templates/order_generic.xlsx')
GLOBAL_SEQUENCING_TEMPLATE = resource_path('sequencing_template',REPO_ROOT/'resources/templates/sequencing_generic.xlsx')
LEGACY_CONFIRMED_SEQUENCING_TEMPLATE = resource_path('legacy_sequencing_template',GLOBAL_SEQUENCING_TEMPLATE)
OPERATIONAL_LABEL_FORMAT = "<minimal_backbone>_<species+canonical_identity_or_species-aware_symbol>"
PRIMER_NAME_FORMAT = "<operational_label>_<F|R>"
SEQUENCING_SAMPLE_FORMAT = "<operational_label>-<colony_number>"

CANONICAL_DIRECTORIES = (
    "inputs",
    "results",
    "results/primers",
    "results/maps",
    "results/order",
    "results/sequencing",
    "results/qc",
    "results/report",
    "scripts",
)

CANONICAL_OUTPUTS = {
    "construct_plan": "construct_plan.tsv",
    "primer_master": "results/primers/primer_master.tsv",
    "primer_qc": "results/qc/primer_qc.tsv",
    "maps": "results/maps",
    "order": "results/order/primer_order.xlsx",
    "sequencing": "results/sequencing/sequencing_order.xlsx",
    "report": "results/report",
}

REQUIRED_PROJECT_FILES = (
    "WORK_LOG.md",
    "project.yaml",
    "inputs/runtime_inputs.json",
    CANONICAL_OUTPUTS["construct_plan"],
    CANONICAL_OUTPUTS["primer_master"],
    CANONICAL_OUTPUTS["primer_qc"],
    CANONICAL_OUTPUTS["order"],
    CANONICAL_OUTPUTS["sequencing"],
)

SPECIES_ABBREVIATIONS = {
    "rice": "Os",
    "oryza sativa": "Os",
    "arabidopsis": "At",
    "arabidopsis thaliana": "At",
    "maize": "Zm",
    "zea mays": "Zm",
    "nicotiana benthamiana": "Nb",
    "tomato": "Sl",
    "solanum lycopersicum": "Sl",
    "marchantia": "Mp",
    "marchantia polymorpha": "Mp",
}

MINIMAL_BACKBONE_ALIASES = {
    "PHW_3xHA": "PHW",
    "PHW": "PHW",
    "pGreenII-62SK": "62SK",
    "62SK": "62SK",
    "pGreenII-0800": "0800",
    "0800": "0800",
    "nYFP": "nYFP",
    "cYFP": "cYFP",
}

PRIMER_MASTER_FIELDS = (
    "project",
    "species",
    "gene_id",
    "symbol",
    "canonical_identity",
    "operational_label",
    "full_construct_name",
    "backbone",
    "preset_id",
    "preset_version",
    "primer_direction",
    "primer_name",
    "annealing_sequence",
    "homology_arm",
    "final_sequence",
    "Tm",
    "GC",
    "warning",
    "selected_reason",
    "remove_stop",
    "physical_plasmid_status",
    "status",
)

CONSTRUCT_PLAN_FIELDS = (
    "project",
    "species",
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
    "identity_decision",
    "overrides",
    "status",
)

PRIMER_QC_FIELDS = (
    "project",
    "gene_id",
    "canonical_identity",
    "full_construct_name",
    "check",
    "status",
    "detail",
    "extra_qc",
)

CANONICAL_QC_CHECKS = (
    "identity_audit",
    "gene_identity",
    "CDS",
    "start_codon",
    "stop_removal",
    "frame",
    "homology_arm",
    "Tm",
    "GC",
    "pair_delta",
    "hairpin",
    "homodimer",
    "heterodimer",
    "fusion",
    "map",
    "order",
    "sequencing",
)

ALLOWED_OVERRIDE_TYPES = {
    "cloning_architecture",
    "fusion",
    "reporter_architecture",
    "sequencing_service",
    "user_requested_name",
    "scientific_qc",
    "additional_report",
    "colony_count",
}

IDENTITY_PREFLIGHT_SOURCE_KINDS = (
    "identity_registry",
    "canonical_construct_master",
    "active_construct_plan",
    "historical_mapping",
    "physical_inventory",
)

# Rice TF/RTF identity assignment is anchored by the current 2026 master.
# Other sources enrich collision and physical-status evidence but must never
# turn an otherwise valid new construction into a global false block merely
# because a redundant historical table does not exist.
IDENTITY_SOURCE_SEVERITY = {
    "identity_registry": "HARD_REQUIRED",
    "canonical_construct_master": "SUPPLEMENTARY",
    "active_construct_plan": "SUPPLEMENTARY",
    "historical_mapping": "SUPPLEMENTARY",
    "physical_inventory": "OPTIONAL_STATUS",
}
HARD_REQUIRED_IDENTITY_SOURCE_KINDS = tuple(
    kind for kind in IDENTITY_PREFLIGHT_SOURCE_KINDS
    if IDENTITY_SOURCE_SEVERITY[kind] == "HARD_REQUIRED"
)


def canonical_project_paths(project_root: str | Path) -> dict[str, Path]:
    """Return every canonical contract path without touching the filesystem."""
    root = Path(project_root)
    paths = {key: root / relative for key, relative in CANONICAL_OUTPUTS.items()}
    paths.update(
        {
            "root": root,
            "work_log": root / "WORK_LOG.md",
            "manifest": root / "project.yaml",
            "inputs": root / "inputs",
            "scripts": root / "scripts",
        }
    )
    return paths


def canonical_species_abbreviation(species: str) -> str:
    """Resolve a supported species name or already-canonical abbreviation."""
    value = str(species or "").strip()
    if value in set(SPECIES_ABBREVIATIONS.values()):
        return value
    abbreviation = SPECIES_ABBREVIATIONS.get(value.lower())
    if not abbreviation:
        raise ValueError(f"Unsupported species without a canonical abbreviation: {species!r}")
    return abbreviation


def minimal_backbone_name(backbone: str) -> str:
    """Resolve the canonical operational backbone display name."""
    value = str(backbone or "").strip()
    alias = MINIMAL_BACKBONE_ALIASES.get(value)
    if alias:return alias
    # An explicit custom short name is a naming label, never an identity/sequence match.
    import re
    if re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',value):return value
    raise ValueError('Provide a simple explicit backbone short name')
