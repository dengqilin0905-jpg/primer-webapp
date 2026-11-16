"""Canonical Rice gene annotation resolution for primer workflows."""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path
from resource_config import resource_path,REPO_ROOT
from typing import Iterable

from naming import normalize_species_symbol


PRIMARY_ANNOTATION_PATH = resource_path('rice_annotation',REPO_ROOT/'resources/examples/annotation.tsv')
CURATED_ANNOTATION_PATH = resource_path('curated_annotation',REPO_ROOT/'resources/examples/curated_annotation.tsv')
PRIMARY_REQUIRED_COLUMNS = {
    "geneID",
    "gene_symbol",
    "AtBestHit",
    "Arabidopsis_symbol",
}
CURATED_REQUIRED_COLUMNS = {"Gene_ID", "Gene_symbol"}
MISSING_VALUES = {"", "-", ".", "na", "n/a", "none", "null", "unknown"}


class RiceAnnotationError(ValueError):
    """Base error for canonical Rice annotation resolution."""


class RiceAnnotationSourceError(RiceAnnotationError):
    """Raised when an annotation source cannot be read."""


class RiceAnnotationSchemaError(RiceAnnotationError):
    """Raised when an annotation source lacks required columns."""


@dataclass(frozen=True)
class RiceAnnotationResult:
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

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def _is_present(value: object) -> bool:
    return str(value or "").strip().lower() not in MISSING_VALUES


def _split_aliases(value: object) -> tuple[str, ...]:
    aliases = []
    for item in str(value or "").replace(",", ";").split(";"):
        candidate = item.strip()
        if _is_present(candidate) and candidate not in aliases:
            aliases.append(candidate)
    return tuple(aliases)


def _read_tsv(
    path: Path,
    *,
    required_columns: set[str],
    source_label: str,
) -> list[dict[str, str]]:
    try:
        handle = path.open(encoding="utf-8-sig", newline="")
    except OSError as exc:
        raise RiceAnnotationSourceError(
            f"Unable to read {source_label} Rice annotation source {path}: {exc}"
        ) from exc
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = set(reader.fieldnames or [])
        missing = sorted(required_columns - columns)
        if missing:
            raise RiceAnnotationSchemaError(
                f"{source_label} Rice annotation source {path} lacks required columns: {missing}"
            )
        return [dict(row) for row in reader]


def _index_rows(
    rows: Iterable[dict[str, str]],
    *,
    gene_column: str,
    source_path: Path,
) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        gene_id = str(row.get(gene_column) or "").strip()
        if not gene_id:
            continue
        if gene_id in indexed and indexed[gene_id] != row:
            raise RiceAnnotationSchemaError(
                f"Conflicting duplicate annotation rows for {gene_id} in {source_path}"
            )
        indexed[gene_id] = row
    return indexed


def _at_homolog(row: dict[str, str] | None) -> tuple[str | None, str | None]:
    if not row:
        return None, None
    symbols = _split_aliases(row.get("Arabidopsis_symbol"))
    at_symbol = symbols[0] if symbols else None
    at_gene_id = None
    best_hit = str(row.get("AtBestHit") or "").strip().split(";", 1)[0].strip()
    if _is_present(best_hit):
        marker = best_hit.rfind("_AT")
        if marker >= 0:
            best_symbols = _split_aliases(best_hit[:marker])
            at_gene_id = best_hit[marker + 1 :]
            if not at_symbol and best_symbols:
                at_symbol = best_symbols[0]
        elif best_hit.upper().startswith("AT"):
            at_gene_id = best_hit
    return at_gene_id, at_symbol


class RiceAnnotationResolver:
    """Resolve Rice display symbols from one authoritative annotation path."""

    def __init__(
        self,
        primary_path: str | Path = PRIMARY_ANNOTATION_PATH,
        *,
        curated_path: str | Path | None = CURATED_ANNOTATION_PATH,
    ) -> None:
        self.primary_path = Path(primary_path)
        primary_rows = _read_tsv(
            self.primary_path,
            required_columns=PRIMARY_REQUIRED_COLUMNS,
            source_label="primary",
        )
        self._primary = _index_rows(
            primary_rows,
            gene_column="geneID",
            source_path=self.primary_path,
        )

        self.curated_path = Path(curated_path) if curated_path is not None else None
        self._curated: dict[str, dict[str, str]] = {}
        if self.curated_path is not None:
            curated_rows = _read_tsv(
                self.curated_path,
                required_columns=CURATED_REQUIRED_COLUMNS,
                source_label="curated",
            )
            self._curated = _index_rows(
                curated_rows,
                gene_column="Gene_ID",
                source_path=self.curated_path,
            )

    def resolve(self, gene_id: str, *, input_hint: str | None = None) -> RiceAnnotationResult:
        stable_gene_id = str(gene_id or "").strip()
        if not stable_gene_id:
            raise RiceAnnotationError("Rice gene ID must not be empty")

        primary = self._primary.get(stable_gene_id)
        curated = self._curated.get(stable_gene_id)
        official_aliases = _split_aliases(primary.get("gene_symbol") if primary else None)
        curated_aliases = _split_aliases(curated.get("Gene_symbol") if curated else None)
        at_gene_id, at_symbol = _at_homolog(primary)
        warnings: list[str] = []

        if official_aliases:
            resolved_symbol = normalize_species_symbol(
                official_aliases[0], species="Rice", official_symbol_available=True
            )
            official_symbol = official_aliases[0]
            source = "PRIMARY_OFFICIAL:gene_symbol"
            source_path = str(self.primary_path)
            fallback_type = "OFFICIAL_RICE_SYMBOL"
            accepted_hints = set(official_aliases)
        elif curated_aliases:
            resolved_symbol = normalize_species_symbol(
                curated_aliases[0], species="Rice", official_symbol_available=True
            )
            official_symbol = curated_aliases[0]
            source = "CURATED:Gene_symbol"
            source_path = str(self.curated_path)
            fallback_type = "CURATED_RICE_SYMBOL"
            accepted_hints = set(curated_aliases)
            warnings.append("PRIMARY_OFFICIAL_SYMBOL_UNAVAILABLE; used curated Rice annotation")
        elif at_symbol:
            resolved_symbol = normalize_species_symbol(
                at_symbol,
                species="Rice",
                at_homolog=at_symbol,
                official_symbol_available=False,
            )
            official_symbol = None
            source = "PRIMARY_HOMOLOG:Arabidopsis_symbol/AtBestHit"
            source_path = str(self.primary_path)
            fallback_type = "ARABIDOPSIS_HOMOLOG"
            accepted_hints = {resolved_symbol, at_symbol}
            warnings.append("OFFICIAL_RICE_SYMBOL_UNAVAILABLE; display inferred from Arabidopsis homolog")
        else:
            resolved_symbol = stable_gene_id
            official_symbol = None
            source = "STABLE_GENE_ID"
            source_path = str(self.primary_path)
            fallback_type = "STABLE_GENE_ID"
            accepted_hints = {stable_gene_id}
            warnings.append("ANNOTATION_UNAVAILABLE; retained stable Rice gene ID")

        raw_hint = str(input_hint or "").strip()
        hint = raw_hint if _is_present(raw_hint) else None
        hint_conflict = bool(hint and hint not in accepted_hints and hint != resolved_symbol)
        if hint_conflict:
            warnings.append(
                f"INPUT_HINT_CONFLICT: {hint} does not match canonical {resolved_symbol}; canonical retained"
            )

        return RiceAnnotationResult(
            gene_id=stable_gene_id,
            resolved_symbol=resolved_symbol,
            official_rice_symbol=official_symbol,
            at_homolog_id=at_gene_id,
            at_homolog_symbol=at_symbol,
            source=source,
            source_path=source_path,
            fallback_type=fallback_type,
            warning="; ".join(warnings),
            input_hint=hint,
            hint_conflict=hint_conflict,
        )


def resolve_rice_annotation(
    gene_id: str,
    *,
    input_hint: str | None = None,
    primary_path: str | Path = PRIMARY_ANNOTATION_PATH,
    curated_path: str | Path | None = CURATED_ANNOTATION_PATH,
) -> RiceAnnotationResult:
    """Convenience resolver for callers that do not need to reuse an index."""
    return RiceAnnotationResolver(
        primary_path,
        curated_path=curated_path,
    ).resolve(gene_id, input_hint=input_hint)
