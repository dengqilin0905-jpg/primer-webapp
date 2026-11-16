"""Canonical business-result assembly for ordinary primer projects.

This module is the only layer that combines annotation, naming, primer design,
lengths, rescue state, exporter inputs, warnings, provenance, and QC rows.
Exporters consume the assembled result and do not resolve those values again.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from typing import Any, Mapping

from design_config import resolve_batch_backbone
from genbank_export import build_genbank_record
from naming import (
    build_construct_name,
    build_map_filename,
    build_operational_label,
    build_primer_label,
)
from primer_design import STOP_CODONS, design_primers_for_gene, normalize_seq
from project_contract import (
    CANONICAL_NORMALIZATION_VERSION,
    CANONICAL_QC_CHECKS,
    CANONICAL_RESULT_VERSION,
    canonical_species_abbreviation,
)
from project_manifest import validate_project_manifest
from rice_annotation import RiceAnnotationResolver
from sequencing_export import DEFAULT_COLONY_COUNT


STANDARD_RUNNER_OVERRIDE_TYPES = {"colony_count", "scientific_qc", "additional_report"}
DESIGN_ALLOWED_IDENTITY_DECISIONS = {
    "REBUILD_EXISTING_IDENTITY",
    "APPROVED_NEW_IDENTITY",
    "NOT_APPLICABLE",
}
STRUCTURAL_QC_CHECKS = ("hairpin", "homodimer", "heterodimer")
QC_STATUSES = {"PASS", "WARNING", "FAIL"}


class ProjectServiceError(ValueError):
    """Raised when canonical result assembly cannot safely continue."""


@dataclass(frozen=True)
class CanonicalAnnotation:
    gene_id: str
    resolved_symbol: str
    official_rice_symbol: str | None
    at_homolog_id: str | None
    at_homolog_symbol: str | None
    source: str
    source_path: str
    fallback_type: str
    warning: str
    input_hint: str | None
    hint_conflict: bool


@dataclass(frozen=True)
class CanonicalPrimerResult:
    direction: str
    name: str
    annealing_sequence: str
    homology_arm: str
    final_sequence: str
    annealing_length: int
    homology_length: int
    final_length: int
    tm: float
    gc: float


@dataclass(frozen=True)
class CanonicalSequenceResult:
    original_cds: str
    insert: str
    backbone: str


@dataclass(frozen=True)
class CanonicalLengthResult:
    original_cds_length: int
    insert_length: int
    amplicon_length: int
    backbone_length: int
    replaced_backbone_length: int
    final_construct_length: int
    sequencing_length: int


@dataclass(frozen=True)
class CanonicalRescueResult:
    rescue_tier: str
    selected_reason: str
    tm_profile: str
    tm_delta: float
    failure_trace: tuple[Mapping[str, Any], ...]


@dataclass(frozen=True)
class CanonicalProvenance:
    annotation_source: str
    annotation_source_path: str
    sequence_source: str
    sequence_sha256: str
    backbone_source: str
    vector_sha256: str
    preset_id: str
    preset_version: int
    primer_engine: str
    naming_service: str
    map_builder: str


@dataclass(frozen=True)
class CanonicalConstructResult:
    project_id: str
    species: str
    gene_id: str
    canonical_identity: str
    annotation: CanonicalAnnotation
    operational_label: str
    full_construct_name: str
    backbone: str
    preset_id: str
    preset_version: int
    remove_stop: bool
    physical_plasmid_status: str
    identity_decision: str
    sequences: CanonicalSequenceResult
    primers: tuple[CanonicalPrimerResult, CanonicalPrimerResult]
    lengths: CanonicalLengthResult
    rescue: CanonicalRescueResult
    map_input: Mapping[str, Any]
    order_input: Mapping[str, Any]
    sequencing_input: Mapping[str, Any]
    primer_master_rows: tuple[Mapping[str, Any], Mapping[str, Any]]
    construct_plan_row: Mapping[str, Any]
    qc_rows: tuple[Mapping[str, str], ...]
    warnings: tuple[str, ...]
    provenance: CanonicalProvenance
    status: str
    colony_count: int

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class CanonicalProjectResult:
    schema_version: int
    project_id: str
    manifest: Mapping[str, Any]
    constructs: tuple[CanonicalConstructResult, ...]
    colony_count: int

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def normalize_canonical_result(result: CanonicalProjectResult) -> dict[str, Any]:
    """Return the sole JSON-stable Web/Batch parity representation.

    The entrypoints must relay this value unchanged.  Business fields are
    selected here so neither Web nor Batch/Codex can independently rebuild
    annotation, names, primer metrics, lengths, exporter inputs, QC, or
    provenance.
    """
    constructs = []
    for item in result.constructs:
        constructs.append(
            {
                "identity": {
                    "gene_id": item.gene_id,
                    "canonical_identity": item.canonical_identity,
                    "identity_decision": item.identity_decision,
                },
                "annotation": asdict(item.annotation),
                "canonical_symbol": item.annotation.resolved_symbol,
                "naming": {
                    "operational_label": item.operational_label,
                    "full_construct_name": item.full_construct_name,
                    "primer_names": {
                        primer.direction: primer.name for primer in item.primers
                    },
                    "map_filename": item.map_input["output_filename"],
                },
                "primers": [asdict(primer) for primer in item.primers],
                "lengths": asdict(item.lengths),
                "rescue": asdict(item.rescue),
                "warnings": list(item.warnings),
                "status": item.status,
                "map_input": dict(item.map_input),
                "order_input": dict(item.order_input),
                "sequencing_input": dict(item.sequencing_input),
                "qc": [dict(row) for row in item.qc_rows],
                "provenance": asdict(item.provenance),
            }
        )
    payload = {
        "normalization_version": CANONICAL_NORMALIZATION_VERSION,
        "canonical_result_version": result.schema_version,
        "project_id": result.project_id,
        "constructs": constructs,
    }
    # JSON round-trip canonicalizes tuples and mapping implementations while
    # rejecting any non-JSON business value instead of silently stringifying it.
    return json.loads(
        json.dumps(
            payload,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
    )


def normalized_canonical_result_sha256(normalized: Mapping[str, Any]) -> str:
    """Hash an already-normalized result without rebuilding any fields."""
    encoded = json.dumps(
        normalized,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _species_for_construct(manifest: Mapping[str, Any], construct: Mapping[str, Any]) -> str:
    species = construct.get("species") or manifest.get("species")
    if isinstance(species, list):
        if len(species) != 1:
            raise ProjectServiceError(
                f"{construct['gene_id']}: multi-species project requires construct-level species"
            )
        species = species[0]
    return str(species)


def _overrides(construct: Mapping[str, Any]) -> dict[str, Any]:
    resolved: dict[str, Any] = {}
    for override in construct.get("overrides", []):
        override_type = override["override_type"]
        if override_type not in STANDARD_RUNNER_OVERRIDE_TYPES:
            raise ProjectServiceError(
                f"{construct['gene_id']}: override {override_type!r} requires a specialized project runner"
            )
        if override_type in resolved:
            raise ProjectServiceError(
                f"{construct['gene_id']}: duplicate override {override_type!r}"
            )
        resolved[override_type] = override["override_value"]
    return resolved


def _identity_preflight(
    construct: Mapping[str, Any], identity_audits: Mapping[str, Mapping[str, Any]]
) -> dict[str, Any]:
    gene_id = str(construct["gene_id"])
    audit = dict(identity_audits.get(gene_id) or {})
    if not audit:
        raise ProjectServiceError(f"{gene_id}: missing identity preflight")
    decision = audit.get("decision")
    if decision == "SKIP_REBUILD":
        raise ProjectServiceError(f"{gene_id}: SKIP_REBUILD; ordinary runner must not rebuild it")
    if decision not in DESIGN_ALLOWED_IDENTITY_DECISIONS:
        raise ProjectServiceError(f"{gene_id}: identity preflight blocks design ({decision})")
    existing = str(audit.get("canonical_identity") or "")
    requested = str(construct.get("canonical_identity") or "")
    if existing and existing != requested:
        from identity_audit import identity_key

        if identity_key(existing) != identity_key(requested):
            raise ProjectServiceError(
                f"{gene_id}: manifest identity {requested!r} conflicts with preflight {existing!r}"
            )
    return audit


def _structural_qc(
    gene_id: str, structural_qc: Mapping[str, Mapping[str, Any]]
) -> dict[str, dict[str, str]]:
    supplied = structural_qc.get(gene_id) or {}
    resolved: dict[str, dict[str, str]] = {}
    for check in STRUCTURAL_QC_CHECKS:
        value = supplied.get(check)
        if not isinstance(value, Mapping):
            raise ProjectServiceError(f"{gene_id}: missing canonical structural QC {check}")
        status = str(value.get("status") or "").upper()
        if status not in QC_STATUSES:
            raise ProjectServiceError(f"{gene_id}: invalid {check} QC status {status!r}")
        resolved[check] = {"status": status, "detail": str(value.get("detail") or "")}
        if status == "FAIL":
            raise ProjectServiceError(f"{gene_id}: {check} QC failed")
    return resolved


def _validate_cds(gene_id: str, sequence: str, remove_stop: bool) -> str:
    normalized = normalize_seq(sequence)
    if not normalized or len(normalized) % 3:
        raise ProjectServiceError(f"{gene_id}: CDS must be non-empty and divisible by 3")
    if not normalized.startswith("ATG"):
        raise ProjectServiceError(f"{gene_id}: CDS must start with ATG")
    if normalized[-3:] not in STOP_CODONS:
        raise ProjectServiceError(f"{gene_id}: CDS must end with a canonical stop codon")
    if remove_stop and (len(normalized) - 3) % 3:
        raise ProjectServiceError(f"{gene_id}: stop-less CDS is out of frame")
    return normalized


def _canonical_annotation(
    species: str,
    construct: Mapping[str, Any],
    rice_resolver: RiceAnnotationResolver | None,
) -> CanonicalAnnotation:
    gene_id = str(construct["gene_id"])
    input_hint = str(construct.get("symbol") or "").strip() or None
    if canonical_species_abbreviation(species) == "Os":
        if rice_resolver is None:
            raise ProjectServiceError("Rice annotation resolver was not initialized")
        resolved = rice_resolver.resolve(gene_id, input_hint=input_hint)
        return CanonicalAnnotation(**resolved.as_dict())
    symbol = input_hint or gene_id
    return CanonicalAnnotation(
        gene_id=gene_id,
        resolved_symbol=symbol,
        official_rice_symbol=None,
        at_homolog_id=None,
        at_homolog_symbol=None,
        source="MANIFEST_NON_RICE:symbol",
        source_path="project.yaml",
        fallback_type="NON_RICE_MANIFEST_SYMBOL",
        warning="",
        input_hint=input_hint,
        hint_conflict=False,
    )


def _tag_from_preset(construct: Mapping[str, Any], config: Mapping[str, Any]) -> str | None:
    explicit = str(construct.get("tag") or "").strip()
    if explicit:
        return explicit
    parts = [part.strip() for part in str(config.get("fusion_order") or "").split("->")]
    parts = [part for part in parts if part and part.lower() != "target"]
    if len(parts) == 1 and parts[0].casefold() != str(construct["backbone"]).casefold():
        return parts[0]
    return None


def _qc_row(
    project_id: str,
    gene_id: str,
    canonical_identity: str,
    full_construct_name: str,
    check: str,
    status: str,
    detail: str,
) -> dict[str, str]:
    return {
        "project": project_id,
        "gene_id": gene_id,
        "canonical_identity": canonical_identity,
        "full_construct_name": full_construct_name,
        "check": check,
        "status": status,
        "detail": detail,
        "extra_qc": "",
    }


def _qc_rows(
    *,
    project_id: str,
    gene_id: str,
    canonical_identity: str,
    full_construct_name: str,
    annotation: CanonicalAnnotation,
    design: Mapping[str, Any],
    identity: Mapping[str, Any],
    structural: Mapping[str, Mapping[str, str]],
    remove_stop: bool,
    operational_label: str,
    colony_count: int,
) -> tuple[Mapping[str, str], ...]:
    sequence = normalize_seq(design["original_seq"])
    insert = normalize_seq(design["insert_seq"])
    gene_status = "WARNING" if annotation.warning else "PASS"
    checks = {
        "identity_audit": ("PASS", str(identity.get("decision"))),
        "gene_identity": (
            gene_status,
            f"gene_id={gene_id}; symbol={annotation.resolved_symbol}; "
            f"annotation_source={annotation.source}; warning={annotation.warning}",
        ),
        "CDS": ("PASS", f"input_length={len(sequence)}; insert_length={len(insert)}"),
        "start_codon": ("PASS", sequence[:3]),
        "stop_removal": ("PASS", f"remove_stop={remove_stop}; insert_terminal={insert[-3:]}"),
        "frame": ("PASS", f"insert_length_mod3={len(insert) % 3}"),
        "homology_arm": ("PASS", "Forward and reverse homology arms match the selected preset"),
        "Tm": (
            "PASS" if 55 <= design["tm_f"] <= 65 and 55 <= design["tm_r"] <= 65 else "WARNING",
            f"F={design['tm_f']:.2f}; R={design['tm_r']:.2f}",
        ),
        "GC": (
            "PASS" if 40 <= design["gc_f"] <= 60 and 40 <= design["gc_r"] <= 60 else "WARNING",
            f"F={design['gc_f']:.2f}; R={design['gc_r']:.2f}",
        ),
        "pair_delta": (
            "PASS" if design["tm_delta"] <= 2 else "WARNING",
            f"Tm_delta={design['tm_delta']:.2f}",
        ),
        "hairpin": (structural["hairpin"]["status"], structural["hairpin"]["detail"]),
        "homodimer": (structural["homodimer"]["status"], structural["homodimer"]["detail"]),
        "heterodimer": (structural["heterodimer"]["status"], structural["heterodimer"]["detail"]),
        "fusion": ("PASS", design.get("frame_check_warnings") or "Preset fusion/frame policy evaluated"),
        "map": ("PASS", f"{full_construct_name}.gbk"),
        "order": ("PASS", f"{operational_label}_F/R"),
        "sequencing": ("PASS", f"{operational_label}-1..{colony_count}"),
    }
    return tuple(
        _qc_row(
            project_id,
            gene_id,
            canonical_identity,
            full_construct_name,
            check,
            *checks[check],
        )
        for check in CANONICAL_QC_CHECKS
    )


def assemble_project_result(
    manifest: Mapping[str, Any],
    *,
    sequences: Mapping[str, str],
    identity_audits: Mapping[str, Mapping[str, Any]],
    structural_qc: Mapping[str, Mapping[str, Any]],
    sequencing_inputs: Mapping[str, Mapping[str, Any]],
    rice_annotation_resolver: RiceAnnotationResolver | None = None,
) -> CanonicalProjectResult:
    """Assemble every downstream business value exactly once."""
    checked = validate_project_manifest(dict(manifest))
    needs_rice = any(
        canonical_species_abbreviation(_species_for_construct(checked, item)) == "Os"
        for item in checked["constructs"]
    )
    resolver = rice_annotation_resolver
    if needs_rice and resolver is None:
        resolver = RiceAnnotationResolver()

    results: list[CanonicalConstructResult] = []
    colony_counts: set[int] = set()
    for construct in checked["constructs"]:
        gene_id = str(construct["gene_id"])
        species = _species_for_construct(checked, construct)
        overrides = _overrides(construct)
        identity = _identity_preflight(construct, identity_audits)
        structural = _structural_qc(gene_id, structural_qc)
        sequence = _validate_cds(gene_id, sequences.get(gene_id, ""), bool(construct["remove_stop"]))
        config = resolve_batch_backbone(
            construct["preset_id"],
            version=construct["preset_version"],
            overrides={"remove_stop_override": bool(construct["remove_stop"])},
        )
        annotation = _canonical_annotation(species, construct, resolver)
        canonical_identity = str(construct.get("canonical_identity") or "")
        operational_label = build_operational_label(
            annotation.resolved_symbol,
            species=species,
            identifier=canonical_identity or None,
            backbone=construct["backbone"],
            # The resolver has already produced the canonical display symbol;
            # naming must consume it without attempting a second annotation fallback.
            official_symbol_available=True,
        )
        full_construct_name = build_construct_name(
            construct["backbone"],
            annotation.resolved_symbol,
            species=species,
            identifier=canonical_identity or None,
            tag=_tag_from_preset(construct, config),
            official_symbol_available=True,
        )
        design = design_primers_for_gene(
            gene_name=full_construct_name,
            seq=sequence,
            homo_f=config["homo_f"],
            homo_r=config["homo_r"],
            remove_stop=config["effective_remove_stop"],
            min_len=18,
            max_len=30,
            tm_min=55,
            tm_max=65,
            target_tm=60,
            gc_min=40,
            gc_max=60,
            Na=50,
            Mg=1.5,
            dNTPs=0.2,
            primer_conc=250,
            seq_type="CDS",
            length_mode="adaptive_soft",
            cds_len_pref=(18, 21),
            tm_priority_weight=2.5,
            tm_pair_delta_max=2,
            tm_delta_mode="soft",
            tm_profile="snapgene_preset",
            cds_recovery_enable=True,
            enable_frame_check=config["primer_engine_frame_check"],
            construct_type=config["construct_type"],
            restriction_site_f=config["restriction_site_f"],
            restriction_site_r=config["restriction_site_r"],
        )
        primer_f = build_primer_label(
            annotation.resolved_symbol,
            "F",
            species=species,
            identifier=canonical_identity or None,
            backbone=construct["backbone"],
        )
        primer_r = build_primer_label(
            annotation.resolved_symbol,
            "R",
            species=species,
            identifier=canonical_identity or None,
            backbone=construct["backbone"],
        )
        warnings = [value for value in (annotation.warning, str(design["warnings"])) if value]
        if construct.get("operational_label") != operational_label:
            warnings.append(
                f"MANIFEST_OPERATIONAL_NAME_DRIFT: {construct.get('operational_label')} -> {operational_label}"
            )
        if construct.get("full_construct_name") != full_construct_name:
            warnings.append(
                f"MANIFEST_CONSTRUCT_NAME_DRIFT: {construct.get('full_construct_name')} -> {full_construct_name}"
            )
        status = "WARNING" if warnings else "PASS"

        primers = (
            CanonicalPrimerResult(
                "F", primer_f, design["gene_bind_f"], config["homo_f"], design["final_f"],
                len(design["gene_bind_f"]), len(config["homo_f"]), len(design["final_f"]),
                float(design["tm_f"]), float(design["gc_f"]),
            ),
            CanonicalPrimerResult(
                "R", primer_r, design["gene_bind_r"], config["homo_r"], design["final_r"],
                len(design["gene_bind_r"]), len(config["homo_r"]), len(design["final_r"]),
                float(design["tm_r"]), float(design["gc_r"]),
            ),
        )
        annotation_note = (
            f"gene_ID={gene_id}; canonical_identity={canonical_identity}; "
            f"canonical_symbol={annotation.resolved_symbol}; official_Os_symbol="
            f"{annotation.official_rice_symbol or 'unavailable'}; At_homolog_ID="
            f"{annotation.at_homolog_id or 'unavailable'}; At_homolog_symbol="
            f"{annotation.at_homolog_symbol or 'unavailable'}; annotation_source="
            f"{annotation.source}; annotation_source_path={annotation.source_path}; "
            f"fallback_type={annotation.fallback_type}"
        )
        map_input = {
            **design,
            "primer_name_f": primer_f,
            "primer_name_r": primer_r,
            "backbone_name": config["backbone_name"],
            "gene_name": annotation.resolved_symbol,
            "plasmid_name": full_construct_name,
            "output_filename": build_map_filename(full_construct_name),
            "homo_f": config["homo_f"],
            "homo_r": config["homo_r"],
            "vector_sequence": config["vector_sequence"],
            "vector_features_json": config["vector_features_json"],
            "vector_topology": config["vector_topology"],
            "vector_source_format": config["vector_source_format"],
            "vector_record_name": config["vector_record_name"],
            "insertion_mode": config["insertion_mode"],
            "restriction_site_f": config["restriction_site_f"],
            "restriction_site_r": config["restriction_site_r"],
            "insert_feature_type": "CDS",
            "insert_feature_name": annotation.resolved_symbol,
            "insert_feature_note": f"Inserted sequence; {annotation_note}",
            "insert_gene": annotation.resolved_symbol,
            "insert_locus_tag": gene_id,
            "vector_map_note": "; ".join(
                value for value in (config.get("vector_map_note", ""), annotation_note) if value
            ),
        }
        record = build_genbank_record(map_input)

        sequencing = dict(sequencing_inputs.get(gene_id) or {})
        direction = str(sequencing.get("reference_primer_direction") or "").upper()
        if direction not in {"F", "R"}:
            raise ProjectServiceError(
                f"{gene_id}: sequencing reference direction must be explicit F/R"
            )
        expected_reference = config["homo_f"] if direction == "F" else config["homo_r"]
        if normalize_seq(sequencing.get("reference_sequence", "")) != expected_reference:
            raise ProjectServiceError(
                f"{gene_id}: sequencing reference sequence does not match preset {direction} arm"
            )
        sequencing_length = min(len(design["insert_seq"]) + len(expected_reference), 1499)
        replacement = config.get("validated_replacement") or {}
        replaced_length = int(replacement.get("boundary_end", 0)) - int(
            replacement.get("boundary_start", 0)
        )
        lengths = CanonicalLengthResult(
            original_cds_length=len(design["original_seq"]),
            insert_length=len(design["insert_seq"]),
            amplicon_length=int(design["amplicon_len"]),
            backbone_length=len(config["vector_sequence"]),
            replaced_backbone_length=max(0, replaced_length),
            final_construct_length=len(record.seq),
            sequencing_length=sequencing_length,
        )
        canonical_sequences = CanonicalSequenceResult(
            original_cds=str(design["original_seq"]),
            insert=str(design["insert_seq"]),
            backbone=str(config["vector_sequence"]),
        )
        order_input = {
            "primer_name_f": primer_f,
            "primer_name_r": primer_r,
            "final_f": design["final_f"],
            "final_r": design["final_r"],
            "warnings": design["warnings"],
            "selected_reason": design["selected_reason"],
            "tm_delta": design["tm_delta"],
            "rescue_tier": design["rescue_tier"],
            "tm_profile": design["tm_profile"],
        }
        sequencing_input = {
            "operational_label": operational_label,
            "full_construct_name": full_construct_name,
            "reference_primer_name": f"{operational_label}_{direction}",
            "reference_sequence": expected_reference,
            "plasmid_length": lengths.final_construct_length,
            "sequencing_length": lengths.sequencing_length,
            "resistance": str(sequencing.get("resistance") or ""),
            "source_map": f"results/maps/{full_construct_name}.gbk",
        }
        primer_rows = tuple(
            {
                "project": checked["project_id"],
                "species": species,
                "gene_id": gene_id,
                "symbol": annotation.resolved_symbol,
                "canonical_identity": canonical_identity,
                "operational_label": operational_label,
                "full_construct_name": full_construct_name,
                "backbone": construct["backbone"],
                "preset_id": construct["preset_id"],
                "preset_version": construct["preset_version"],
                "primer_direction": primer.direction,
                "primer_name": primer.name,
                "annealing_sequence": primer.annealing_sequence,
                "homology_arm": primer.homology_arm,
                "final_sequence": primer.final_sequence,
                "Tm": primer.tm,
                "GC": primer.gc,
                "warning": design["warnings"],
                "selected_reason": design["selected_reason"],
                "remove_stop": construct["remove_stop"],
                "physical_plasmid_status": construct["physical_plasmid_status"],
                "status": status,
            }
            for primer in primers
        )
        plan_row = {
            "project": checked["project_id"],
            "species": species,
            "gene_id": gene_id,
            "symbol": annotation.resolved_symbol,
            "canonical_identity": canonical_identity,
            "operational_label": operational_label,
            "full_construct_name": full_construct_name,
            "backbone": construct["backbone"],
            "preset_id": construct["preset_id"],
            "preset_version": construct["preset_version"],
            "remove_stop": construct["remove_stop"],
            "physical_plasmid_status": construct["physical_plasmid_status"],
            "identity_decision": identity["decision"],
            "overrides": json.dumps(construct.get("overrides", []), ensure_ascii=False, sort_keys=True),
            "status": status,
        }
        colony_count = int(overrides.get("colony_count", DEFAULT_COLONY_COUNT))
        colony_counts.add(colony_count)
        qc_rows = _qc_rows(
            project_id=checked["project_id"],
            gene_id=gene_id,
            canonical_identity=canonical_identity,
            full_construct_name=full_construct_name,
            annotation=annotation,
            design=design,
            identity=identity,
            structural=structural,
            remove_stop=bool(construct["remove_stop"]),
            operational_label=operational_label,
            colony_count=colony_count,
        )
        provenance = CanonicalProvenance(
            annotation_source=annotation.source,
            annotation_source_path=annotation.source_path,
            sequence_source="runtime_inputs.sequences",
            sequence_sha256=hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            backbone_source=str(config["backbone_source"]),
            vector_sha256=str(config["vector_sha256"]),
            preset_id=str(config["preset_id"]),
            preset_version=int(config["preset_version"]),
            primer_engine="primer_design.design_primers_for_gene",
            naming_service="naming.build_construct_name/build_operational_label/build_primer_label",
            map_builder="genbank_export.build_genbank_record",
        )
        rescue = CanonicalRescueResult(
            rescue_tier=str(design["rescue_tier"]),
            selected_reason=str(design["selected_reason"]),
            tm_profile=str(design["tm_profile"]),
            tm_delta=float(design["tm_delta"]),
            failure_trace=tuple(design["failure_trace"]),
        )
        results.append(
            CanonicalConstructResult(
                project_id=checked["project_id"],
                species=species,
                gene_id=gene_id,
                canonical_identity=canonical_identity,
                annotation=annotation,
                operational_label=operational_label,
                full_construct_name=full_construct_name,
                backbone=str(construct["backbone"]),
                preset_id=str(construct["preset_id"]),
                preset_version=int(construct["preset_version"]),
                remove_stop=bool(construct["remove_stop"]),
                physical_plasmid_status=str(construct["physical_plasmid_status"]),
                identity_decision=str(identity["decision"]),
                sequences=canonical_sequences,
                primers=primers,
                lengths=lengths,
                rescue=rescue,
                map_input=map_input,
                order_input=order_input,
                sequencing_input=sequencing_input,
                primer_master_rows=primer_rows,
                construct_plan_row=plan_row,
                qc_rows=qc_rows,
                warnings=tuple(warnings),
                provenance=provenance,
                status=status,
                colony_count=colony_count,
            )
        )

    if len(colony_counts) != 1:
        raise ProjectServiceError("Ordinary project requires one project-wide colony_count")
    return CanonicalProjectResult(
        schema_version=CANONICAL_RESULT_VERSION,
        project_id=checked["project_id"],
        manifest=checked,
        constructs=tuple(results),
        colony_count=colony_counts.pop(),
    )
