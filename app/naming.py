"""Shared, explicit naming rules for new primer constructions."""

from __future__ import annotations

import re
import unicodedata

from project_contract import SPECIES_ABBREVIATIONS, canonical_species_abbreviation, minimal_backbone_name


def _safe_component(value: str | None, fallback: str = "unknown") -> str:
    text = unicodedata.normalize("NFKC", str(value or "").strip())
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^A-Za-z0-9_.-]", "_", text)
    text = re.sub(r"_+", "_", text).strip("_.")
    return text or fallback


def normalize_species_symbol(
    symbol: str | None = None,
    *,
    species: str | None = None,
    at_homolog: str | None = None,
    official_symbol_available: bool | None = None,
) -> str:
    """Return an explicit display symbol without inventing an Os symbol.

    For rice, callers must set ``official_symbol_available=False`` when the
    supplied symbol is only an Arabidopsis homolog. This makes the inferred
    ``Os_At_<AtSymbol>`` status explicit rather than guessed from free text.
    """
    species_text = str(species or "").strip().lower()
    is_rice = species_text in {"os", "rice", "oryza", "oryza sativa"} or species_text.startswith("oryza sativa")
    if is_rice and official_symbol_available is False:
        homolog = _safe_component(at_homolog or symbol)
        if not homolog:
            raise ValueError("Rice inferred naming requires an Arabidopsis homolog symbol")
        return f"Os_At_{homolog}"
    display = _safe_component(symbol)
    if is_rice and official_symbol_available is None and not display.startswith("Os") and at_homolog:
        raise ValueError("Set official_symbol_available explicitly for a rice symbol")
    return display


def build_construct_name(
    backbone: str,
    symbol: str | None = None,
    *,
    species: str | None = None,
    at_homolog: str | None = None,
    official_symbol_available: bool | None = None,
    identifier: str | None = None,
    tag: str | None = None,
) -> str:
    """Build ``backbone[_identifier]_<display_symbol>[_tag]`` explicitly."""
    parts = [_safe_component(backbone)]
    if identifier:
        parts.append(_safe_component(identifier))
    parts.append(
        normalize_species_symbol(
            symbol,
            species=species,
            at_homolog=at_homolog,
            official_symbol_available=official_symbol_available,
        )
    )
    if tag:
        parts.append(_safe_component(tag))
    return "_".join(part for part in parts if part)


def build_primer_name(construct_name: str, direction: str) -> str:
    """Return the legacy full-construct primer name for historical compatibility."""
    suffix = "F" if str(direction).upper().startswith("F") else "R"
    return f"{_safe_component(construct_name)}_{suffix}"


def _species_prefix_from_symbol(symbol: str | None) -> str | None:
    value = _safe_component(symbol)
    for abbreviation in sorted(set(SPECIES_ABBREVIATIONS.values())):
        remainder = value[len(abbreviation):]
        if value.startswith(abbreviation) and remainder and (remainder[0].isupper() or remainder[0].isdigit() or remainder[0] == "_"):
            return abbreviation
    return None


def _species_aware_identity(identity: str, species_abbreviation: str) -> str:
    value = _safe_component(identity)
    existing_prefix = _species_prefix_from_symbol(value)
    if existing_prefix and existing_prefix != species_abbreviation:
        raise ValueError(
            f"Canonical identity species prefix {existing_prefix!r} conflicts with {species_abbreviation!r}"
        )
    return value if existing_prefix else f"{species_abbreviation}{value}"


def build_operational_label(
    symbol: str | None,
    *,
    species: str | None = None,
    identifier: str | None = None,
    backbone: str,
    at_homolog: str | None = None,
    official_symbol_available: bool | None = None,
) -> str:
    """Build the canonical operational construct label.

    The label is strictly ``minimal backbone + species + canonical identity``
    or, when no identity exists, ``minimal backbone + species-aware symbol``.
    It intentionally excludes gene IDs, tags, and other construct metadata.
    """
    inferred_species = _species_prefix_from_symbol(symbol)
    if species is None:
        if not inferred_species:
            raise ValueError("Operational naming requires species or a species-aware symbol")
        species_abbreviation = inferred_species
    else:
        species_abbreviation = canonical_species_abbreviation(species)
        if inferred_species and inferred_species != species_abbreviation:
            raise ValueError(
                f"Symbol species prefix {inferred_species!r} conflicts with {species_abbreviation!r}"
            )

    if identifier:
        identity_or_symbol = _species_aware_identity(identifier, species_abbreviation)
    else:
        display_symbol = normalize_species_symbol(
            symbol,
            species=species_abbreviation,
            at_homolog=at_homolog,
            official_symbol_available=official_symbol_available,
        )
        symbol_prefix = _species_prefix_from_symbol(display_symbol)
        if symbol_prefix and symbol_prefix != species_abbreviation:
            raise ValueError(
                f"Display symbol species prefix {symbol_prefix!r} conflicts with {species_abbreviation!r}"
            )
        identity_or_symbol = display_symbol if symbol_prefix else f"{species_abbreviation}{display_symbol}"

    return f"{minimal_backbone_name(backbone)}_{identity_or_symbol}"


def build_primer_label(
    symbol: str | None,
    direction: str,
    *,
    species: str | None = None,
    identifier: str | None = None,
    backbone: str,
    at_homolog: str | None = None,
    official_symbol_available: bool | None = None,
) -> str:
    """Return ``<operational_label>_<F|R>`` for primer/order use."""
    suffix = str(direction or "").strip().upper()
    if suffix not in {"F", "R"}:
        raise ValueError(f"Primer direction must be F or R, got {direction!r}")
    operational_label = build_operational_label(
        symbol,
        species=species,
        identifier=identifier,
        backbone=backbone,
        at_homolog=at_homolog,
        official_symbol_available=official_symbol_available,
    )
    return f"{operational_label}_{suffix}"


def build_sequencing_sample_name(operational_label: str, colony_number: int) -> str:
    """Return the canonical fourth-generation sequencing sample name."""
    label = _safe_component(operational_label)
    if label != str(operational_label or "").strip():
        raise ValueError(f"Invalid operational label: {operational_label!r}")
    if not isinstance(colony_number, int) or isinstance(colony_number, bool) or colony_number < 1:
        raise ValueError("colony_number must be a positive integer")
    return f"{label}-{colony_number}"


def build_map_filename(construct_name: str, extension: str = ".gbk") -> str:
    ext = extension if str(extension).startswith(".") else f".{extension}"
    return f"{_safe_component(construct_name)}{ext}"


def construct_from_primer_name(primer_name: str) -> str:
    """Recover a historical construct prefix for map re-export only."""
    value = str(primer_name or "").strip()
    if len(value) > 2 and value[-2] in "_-" and value[-1].upper() in {"F", "R"}:
        return value[:-2]
    return value
