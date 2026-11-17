"""Shared fourth-generation sequencing workbook resolver and exporter."""

from __future__ import annotations

import hashlib
import io
import re
from copy import copy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from openpyxl import load_workbook

from naming import build_sequencing_sample_name
from project_contract import GLOBAL_SEQUENCING_TEMPLATE, LEGACY_CONFIRMED_SEQUENCING_TEMPLATE


UPLOAD_SHEET = "上传模版"
PROTECTED_SHEET = "请勿删除！"
START_ROW = 4
EXPECTED_HEADERS = (
    "* 样品名称",
    "* 样品类型",
    "* 测序项目",
    "* 参考序列或引物序列",
    "* 质粒总长度",
    "* 片段长度",
    "* 抗性",
    "是否返样",
    "返样类型",
    "备注",
)
DEFAULT_COLONY_COUNT = 3
DEFAULT_SAMPLE_TYPE = "平板"
DEFAULT_SERVICE = "插入片段测通"
DEFAULT_RETURN_FLAG = "是"
DEFAULT_RETURN_TYPE = "质粒+菌液"


class SequencingExportError(ValueError):
    """Raised when template, input, or generated workbook violates contract."""


@dataclass(frozen=True)
class SequencingTemplateResolution:
    template_bytes: bytes
    source_path: str
    sha256: str
    source_kind: str


def _workbook_from_bytes(data: bytes):
    try:
        return load_workbook(io.BytesIO(data), data_only=False)
    except Exception as exc:
        raise SequencingExportError(f"Unable to read sequencing workbook: {exc}") from exc


def _validate_template_structure(data: bytes) -> dict[str, Any]:
    workbook = _workbook_from_bytes(data)
    if UPLOAD_SHEET not in workbook.sheetnames:
        raise SequencingExportError(f"Template is missing sheet {UPLOAD_SHEET!r}")
    if PROTECTED_SHEET not in workbook.sheetnames:
        raise SequencingExportError(f"Template is missing sheet {PROTECTED_SHEET!r}")
    if workbook[PROTECTED_SHEET].sheet_state != "hidden":
        raise SequencingExportError(f"Template sheet {PROTECTED_SHEET!r} must remain hidden")
    sheet = workbook[UPLOAD_SHEET]
    headers = tuple(sheet.cell(3, column).value for column in range(1, 11))
    if headers != EXPECTED_HEADERS:
        raise SequencingExportError(f"Unexpected sequencing template headers: {headers!r}")
    return {
        "sheet_names": list(workbook.sheetnames),
        "protected_sheet_state": workbook[PROTECTED_SHEET].sheet_state,
        "validation_count": len(sheet.data_validations.dataValidation),
    }


def resolve_sequencing_template(
    explicit_path: str | Path | None = None,
) -> SequencingTemplateResolution:
    """Resolve one registered template without searching historical projects."""
    if explicit_path is not None:
        path = Path(explicit_path)
        source_kind = "explicit"
    elif GLOBAL_SEQUENCING_TEMPLATE.is_file():
        path = GLOBAL_SEQUENCING_TEMPLATE
        source_kind = "global"
    else:
        path = LEGACY_CONFIRMED_SEQUENCING_TEMPLATE
        source_kind = "registered_legacy_fallback"
    if not path.is_file():
        raise SequencingExportError(f"Sequencing template does not exist: {path}")
    data = path.read_bytes()
    _validate_template_structure(data)
    return SequencingTemplateResolution(
        template_bytes=data,
        source_path=str(path.resolve()),
        sha256=hashlib.sha256(data).hexdigest(),
        source_kind=source_kind,
    )


def _positive_int(value: Any, field: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 1:
        raise SequencingExportError(f"{field} must be a positive integer")
    return value


def _reference_sequence(value: Any) -> str:
    sequence = re.sub(r"\s+", "", str(value or "")).upper()
    if not sequence or set(sequence) - set("ACGT"):
        raise SequencingExportError("reference_sequence must be non-empty A/C/G/T DNA")
    return sequence


def validate_sequencing_items(items: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    """Validate explicit sequence provenance without deriving sequence from names."""
    cleaned: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, raw in enumerate(items, start=1):
        item = dict(raw)
        label = str(item.get("operational_label") or "").strip()
        if not label:
            raise SequencingExportError(f"Item {index}: missing operational_label")
        if label in seen:
            raise SequencingExportError(f"Item {index}: duplicate operational_label {label!r}")
        seen.add(label)
        reference_primer_name = str(item.get("reference_primer_name") or "").strip()
        if reference_primer_name not in {f"{label}_F", f"{label}_R"}:
            raise SequencingExportError(
                f"Item {index}: reference_primer_name must be {label}_F or {label}_R"
            )
        cleaned.append(
            {
                "operational_label": label,
                "full_construct_name": str(item.get("full_construct_name") or "").strip(),
                "reference_primer_name": reference_primer_name,
                "reference_sequence": _reference_sequence(item.get("reference_sequence")),
                "plasmid_length": _positive_int(item.get("plasmid_length"), "plasmid_length"),
                "sequencing_length": _positive_int(item.get("sequencing_length"), "sequencing_length"),
                "resistance": str(item.get("resistance") or "").strip(),
                "source_map": str(item.get("source_map") or "").strip(),
                "remarks": str(item.get("remarks") or "").strip(),
            }
        )
        if not cleaned[-1]["resistance"]:
            raise SequencingExportError(f"Item {index}: missing resistance")
    if not cleaned:
        raise SequencingExportError("At least one sequencing item is required")
    return cleaned


def _copy_row_style(sheet, source_row: int, target_row: int) -> None:
    if sheet.row_dimensions[source_row].height:
        sheet.row_dimensions[target_row].height = sheet.row_dimensions[source_row].height
    for column in range(1, 11):
        source = sheet.cell(source_row, column)
        target = sheet.cell(target_row, column)
        if source.has_style:
            target._style = copy(source._style)
        target.number_format = source.number_format
        target.font = copy(source.font)
        target.fill = copy(source.fill)
        target.border = copy(source.border)
        target.alignment = copy(source.alignment)
        target.protection = copy(source.protection)



def antibiotic_options(workbook):
    definition=workbook.defined_names.get('Antibiotic')
    if definition is None:raise SequencingExportError('Template lacks Antibiotic named range')
    options=[]
    for sheet,area in definition.destinations:
        if sheet not in workbook:raise SequencingExportError('Invalid Antibiotic destination')
        for row in workbook[sheet][area]:
            for cell in row:
                if cell.value is not None and str(cell.value).strip():options.append(str(cell.value).strip())
    if not options:raise SequencingExportError('Antibiotic named range is empty')
    return tuple(dict.fromkeys(options))


def normalize_resistance(value,options):
    value=str(value or '').strip()
    if value in options:return value
    aliases={'ampr':'氨苄','ampicillin':'氨苄','kanr':'卡那','kanamycin':'卡那'}
    canonical=aliases.get(value.casefold())
    if canonical in options:return canonical
    raise SequencingExportError(f'Invalid resistance {value!r}; select a value from template Antibiotic')


def _extend_antibiotic_validation(sheet,last_row):
    from openpyxl.worksheet.cell_range import CellRange
    for dv in sheet.data_validations.dataValidation:
        if str(dv.formula1).lstrip('=')!='Antibiotic':continue
        ranges=[]
        for area in dv.ranges.ranges:
            area=CellRange(str(area))
            if area.min_col==7 and area.max_col==7 and area.max_row>area.min_row:
                area.max_row=max(area.max_row,last_row)
            ranges.append(str(area))
        dv.sqref=' '.join(ranges)

def export_sequencing_workbook(
    items: Iterable[dict[str, Any]],
    *,
    template_bytes: bytes | None = None,
    explicit_template_path: str | Path | None = None,
    colony_count: int = DEFAULT_COLONY_COUNT,
) -> tuple[bytes, dict[str, Any]]:
    """Export canonical operational samples while preserving supplied sequences."""
    colonies = _positive_int(colony_count, "colony_count")
    cleaned = validate_sequencing_items(items)
    resolution = None
    if template_bytes is None:
        resolution = resolve_sequencing_template(explicit_template_path)
        template_bytes = resolution.template_bytes
    template_metadata = _validate_template_structure(template_bytes)
    workbook = _workbook_from_bytes(template_bytes)
    options=antibiotic_options(workbook)
    for item in cleaned:item['resistance']=normalize_resistance(item['resistance'],options)
    sheet = workbook[UPLOAD_SHEET]
    existing_max_row = sheet.max_row
    for row in range(START_ROW, existing_max_row + 1):
        for column in range(1, 11):
            sheet.cell(row, column).value = None

    total_rows = len(cleaned) * colonies
    needed_last_row = START_ROW + total_rows - 1
    for row in range(existing_max_row + 1, needed_last_row + 1):
        _copy_row_style(sheet, START_ROW, row)

    _extend_antibiotic_validation(sheet,needed_last_row)
    row_number = START_ROW
    sample_names: list[str] = []
    for item in cleaned:
        for colony_number in range(1, colonies + 1):
            sample_name = build_sequencing_sample_name(item["operational_label"], colony_number)
            values = (
                sample_name,
                DEFAULT_SAMPLE_TYPE,
                DEFAULT_SERVICE,
                item["reference_sequence"],
                item["plasmid_length"],
                item["sequencing_length"],
                item["resistance"],
                DEFAULT_RETURN_FLAG,
                DEFAULT_RETURN_TYPE,
                item["remarks"] or None,
            )
            for column, value in enumerate(values, start=1):
                sheet.cell(row_number, column).value = value
            sample_names.append(sample_name)
            row_number += 1

    output = io.BytesIO()
    workbook.save(output)
    data = output.getvalue()
    validation = validate_sequencing_workbook(data, cleaned, colony_count=colonies)
    return data, {
        **validation,
        **template_metadata,
        "template_source": resolution.source_path if resolution else "provided_bytes",
        "template_sha256": resolution.sha256 if resolution else hashlib.sha256(template_bytes).hexdigest(),
        "sample_names": sample_names,
        "reference_primers": [item["reference_primer_name"] for item in cleaned],
    }


def validate_sequencing_workbook(
    workbook_bytes: bytes,
    expected_items: Iterable[dict[str, Any]],
    *,
    colony_count: int = DEFAULT_COLONY_COUNT,
) -> dict[str, Any]:
    """Validate all generated rows and template invariants."""
    colonies = _positive_int(colony_count, "colony_count")
    expected = validate_sequencing_items(expected_items)
    structure = _validate_template_structure(workbook_bytes)
    workbook = _workbook_from_bytes(workbook_bytes)
    options=antibiotic_options(workbook)
    for item in expected:item['resistance']=normalize_resistance(item['resistance'],options)
    sheet = workbook[UPLOAD_SHEET]
    observed = []
    for row in range(START_ROW, sheet.max_row + 1):
        values = [sheet.cell(row, column).value for column in range(1, 11)]
        if any(value is not None for value in values):
            if values[6] not in options:raise SequencingExportError('Exported resistance is outside Antibiotic named range')
            coordinate=sheet.cell(row,7).coordinate
            if not any(str(dv.formula1).lstrip('=')=='Antibiotic' and coordinate in dv for dv in sheet.data_validations.dataValidation):
                raise SequencingExportError(f'Missing Antibiotic validation for {coordinate}')
            observed.append(values)
    if len(observed) != len(expected) * colonies:
        raise SequencingExportError(
            f"Sequencing row count {len(observed)} != expected {len(expected) * colonies}"
        )
    by_label = {item["operational_label"]: item for item in expected}
    expected_samples = {
        build_sequencing_sample_name(label, colony)
        for label in by_label
        for colony in range(1, colonies + 1)
    }
    seen = set()
    for values in observed:
        sample_name = values[0]
        if sample_name not in expected_samples or sample_name in seen:
            raise SequencingExportError(f"Unexpected or duplicate sequencing sample: {sample_name!r}")
        seen.add(sample_name)
        label = sample_name.rsplit("-", 1)[0]
        item = by_label[label]
        expected_values = [
            DEFAULT_SAMPLE_TYPE,
            DEFAULT_SERVICE,
            item["reference_sequence"],
            item["plasmid_length"],
            item["sequencing_length"],
            item["resistance"],
            DEFAULT_RETURN_FLAG,
            DEFAULT_RETURN_TYPE,
        ]
        if values[1:9] != expected_values:
            raise SequencingExportError(f"Sequencing workbook values differ for {sample_name}")
    if seen != expected_samples:
        raise SequencingExportError("Sequencing sample set is incomplete")
    return {
        **structure,
        "status": "PASS",
        "construct_count": len(expected),
        "colony_count": colonies,
        "sample_count": len(observed),
    }
