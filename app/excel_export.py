import io
import re
from copy import copy
from pathlib import Path
import openpyxl
from openpyxl.utils import get_column_letter

from resource_config import resource_path,REPO_ROOT
TEMPLATE_PATH = resource_path('order_template',REPO_ROOT/'resources/templates/order_generic.xlsx')
START_ROW = 18
BASE_TEMPLATE_MAX_ROW = 223
# Shared naming.py emits the canonical ``<construct>_F/R`` suffix.  Keep
# compatibility with historical ``<construct>-F/R`` names accepted by older
# exports; both forms must pass validation for current and legacy runs.
PRIMER_NAME_F_RE = re.compile(r"^[A-Za-z0-9_.-]+(?:_|-)F$")
PRIMER_NAME_R_RE = re.compile(r"^[A-Za-z0-9_.-]+(?:_|-)R$")


def _len_formula(row: int) -> str:
    return f"=LEN(C{row})"


def _get_main_sheet(wb):
    if "寡核苷酸合成订购单" in wb.sheetnames:
        return wb["寡核苷酸合成订购单"]
    return wb[wb.sheetnames[0]]


HEADER_DEFAULT_VALUE = 0
HEADER_PROFILE_CELLS = {
        "orderer_name": "B4",
        "mobile": "B5",
        "email": "B6",
        "legal_person": "B7",
        "phone": "B8",
        "institution": "B9",
        "invoice_title": "B10",
        "payment_method": "B11",
        "invoice_with_goods": "B12",
        "shipping_address": "B13",
        "header_note": "B14",
}


def _cell_or_merged_anchor(ws, coordinate: str):
    """Return a writable cell without changing merged-cell layout."""
    for merged_range in ws.merged_cells.ranges:
        if coordinate in merged_range:
            return ws.cell(row=merged_range.min_row, column=merged_range.min_col)
    return ws[coordinate]


def _set_cell_or_merged_anchor(ws, coordinate: str, value) -> str:
    """Write to a normal cell or the top-left anchor of its merged range."""
    cell = _cell_or_merged_anchor(ws, coordinate)
    cell.value = value
    return cell.coordinate


def _clear_cell_or_merged_anchor(ws, coordinate: str) -> str:
    """Clear a normal cell or merged-range anchor while preserving layout."""
    cell = _cell_or_merged_anchor(ws, coordinate)
    cell.value = None
    cell.hyperlink = None
    cell.comment = None
    return cell.coordinate


def _fill_header(ws, profile: dict | None, blank_profile_fields: bool = False):
    default_value = None if blank_profile_fields else HEADER_DEFAULT_VALUE
    for cell in HEADER_PROFILE_CELLS.values():
        _set_cell_or_merged_anchor(ws, cell, default_value)
        if blank_profile_fields:
            # Some uploaded Sangon templates encode contact text as a hyperlink
            # display value; remove the relationship as well as the cell value.
            _clear_cell_or_merged_anchor(ws, cell)

    for k, cell in HEADER_PROFILE_CELLS.items():
        v = profile.get(k) if profile else None
        if v is not None and v != "":
            _set_cell_or_merged_anchor(ws, cell, v)


def _header_text(value) -> str:
    return re.sub(r"\s+", "", str(value or "")).lower()


def _resolve_order_columns(ws) -> dict[str, str]:
    """Resolve order columns from row-17 labels with legacy fallbacks."""
    labels = {
        get_column_letter(column): _header_text(ws.cell(row=17, column=column).value)
        for column in range(1, max(13, ws.max_column) + 1)
    }

    def find(*parts: str) -> str | None:
        for column, label in labels.items():
            if all(part in label for part in parts):
                return column
        return None

    columns = {
        "id": find("id") or "A",
        "name": find("primer", "名称") or "B",
        "sequence": find("序列") or "C",
        "length": find("碱基数") or "D",
        "tubes": find("分装管") or "E",
        "amount": find("提供总量") or "F",
        "unit": find("计量单位") or None,
        "purification": find("纯化") or None,
        "scale": find("nmoles") or None,
        "mod5": find("5'端修饰") or "I",
        "mod3": find("3'端修饰") or "J",
        "dual_mod": find("双标记修饰") or "K",
        "type": find("类型") or "L",
        "remarks": find("备注") or "M",
    }

    # The old built-in workbook has H=nmoles and G=purification. The current
    # standard workbook has F=amount, G=unit, and H=purification.
    columns["purification"] = columns["purification"] or "G"
    columns["scale"] = columns["scale"] or None
    return columns


def _set_delivery_format(ws, delivery_format: str) -> str:
    """Set delivery format through K16, resolving its merged-cell anchor."""
    return _set_cell_or_merged_anchor(ws, "K16", str(delivery_format))


def _copy_row_style(ws, src_row: int, dst_row: int, max_col: int = 13):
    if ws.row_dimensions.get(src_row) and ws.row_dimensions[src_row].height:
        ws.row_dimensions[dst_row].height = ws.row_dimensions[src_row].height
    for c in range(1, max_col + 1):
        src = ws.cell(row=src_row, column=c)
        dst = ws.cell(row=dst_row, column=c)
        if src.has_style:
            dst._style = copy(src._style)
        dst.number_format = src.number_format
        dst.font = copy(src.font)
        dst.fill = copy(src.fill)
        dst.border = copy(src.border)
        dst.alignment = copy(src.alignment)
        dst.protection = copy(src.protection)
        dst.comment = None


def _extend_validations(ws, new_max_row: int):
    try:
        dvs = list(ws.data_validations.dataValidation)
    except Exception:
        return
    if not dvs:
        return

    def _extend_sqref(s: str) -> str:
        parts = str(s).split()
        new_parts = []
        for p in parts:
            m = re.match(r"^([A-Z]+)(\d+):([A-Z]+)(\d+)$", p)
            if m:
                c1, r1, c2, r2 = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
                if r1 == START_ROW and r2 == BASE_TEMPLATE_MAX_ROW and c1 == c2:
                    new_parts.append(f"{c1}{START_ROW}:{c2}{new_max_row}")
                else:
                    new_parts.append(p)
            else:
                new_parts.append(p)
        return " ".join(new_parts)

    for dv in dvs:
        try:
            dv.sqref = _extend_sqref(dv.sqref)
        except Exception:
            pass


def _ensure_rows(ws, needed_max_row: int):
    if needed_max_row <= BASE_TEMPLATE_MAX_ROW:
        return
    _extend_validations(ws, needed_max_row)
    src = BASE_TEMPLATE_MAX_ROW
    for row in range(BASE_TEMPLATE_MAX_ROW + 1, needed_max_row + 1):
        _copy_row_style(ws, src, row)
        ws[f"D{row}"].value = _len_formula(row)


def validate_primer_rows(primer_rows):
    errors = []
    cleaned = []
    for i, r in enumerate(primer_rows, start=1):
        pf = (r.get("primer_name_f") or "").strip()
        pr = (r.get("primer_name_r") or "").strip()
        sf = re.sub(r"[^ACGT]", "", (r.get("final_f") or "").upper())
        sr = re.sub(r"[^ACGT]", "", (r.get("final_r") or "").upper())
        warn = (r.get("warnings") or "").strip()
        selected_reason = (r.get("selected_reason") or "").strip()
        tm_delta = r.get("tm_delta")
        rescue_tier = (r.get("rescue_tier") or "").strip()
        tm_profile = (r.get("tm_profile") or "").strip()

        if not pf or not pr:
            errors.append(f"Row {i}: missing primer name F/R")
        if not sf or not sr:
            errors.append(f"Row {i}: missing final sequence F/R")
        if pf and not PRIMER_NAME_F_RE.fullmatch(pf):
            errors.append(f"Row {i}: invalid primer_name_f format ({pf})")
        if pr and not PRIMER_NAME_R_RE.fullmatch(pr):
            errors.append(f"Row {i}: invalid primer_name_r format ({pr})")

        cleaned.append(
            {
                "primer_name_f": pf,
                "primer_name_r": pr,
                "final_f": sf,
                "final_r": sr,
                "warnings": warn,
                "selected_reason": selected_reason,
                "tm_delta": tm_delta,
                "rescue_tier": rescue_tier,
                "tm_profile": tm_profile,
            }
        )

    if errors:
        raise ValueError("; ".join(errors[:8]))
    return cleaned


def _export_with_workbook(
    wb,
    primer_rows,
    purification: str = "HAP",
    nmoles: int = 5,
    tubes: int = 1,
    profile: dict | None = None,
    delivery_format: str | None = None,
    blank_profile_fields: bool = False,
):
    cleaned_rows = validate_primer_rows(primer_rows)
    ws = _get_main_sheet(wb)

    _fill_header(ws, profile, blank_profile_fields=blank_profile_fields)
    if delivery_format is not None:
        _set_delivery_format(ws, delivery_format)

    columns = _resolve_order_columns(ws)

    flat = []
    for r in cleaned_rows:
        warn = r["warnings"]
        if r.get("selected_reason"):
            warn = f"{warn} | reason:{r['selected_reason']}" if warn else f"reason:{r['selected_reason']}"
        if r.get("tm_delta") is not None and str(r.get("tm_delta")) != "":
            try:
                dtm = float(r["tm_delta"])
                warn = f"{warn} | dTm:{dtm:.2f}" if warn else f"dTm:{dtm:.2f}"
            except Exception:
                pass
        if r.get("rescue_tier"):
            warn = f"{warn} | {r['rescue_tier']}" if warn else str(r["rescue_tier"])
        if r.get("tm_profile"):
            warn = f"{warn} | profile:{r['tm_profile']}" if warn else f"profile:{r['tm_profile']}"
        flat.append((r["primer_name_f"], r["final_f"], warn))
        flat.append((r["primer_name_r"], r["final_r"], warn))

    needed_max_row = START_ROW + len(flat) - 1 if flat else START_ROW
    _ensure_rows(ws, needed_max_row)

    clear_max_row = max(BASE_TEMPLATE_MAX_ROW, needed_max_row, ws.max_row)
    for row in range(START_ROW, clear_max_row + 1):
        _set_cell_or_merged_anchor(ws, f"{columns['id']}{row}", None)
        _set_cell_or_merged_anchor(ws, f"{columns['name']}{row}", None)
        _set_cell_or_merged_anchor(ws, f"{columns['sequence']}{row}", None)
        _set_cell_or_merged_anchor(ws, f"{columns['length']}{row}", _len_formula(row))
        _set_cell_or_merged_anchor(ws, f"{columns['tubes']}{row}", None)
        _set_cell_or_merged_anchor(ws, f"{columns['amount']}{row}", None)
        if columns["unit"]:
            _set_cell_or_merged_anchor(ws, f"{columns['unit']}{row}", None)
        _set_cell_or_merged_anchor(ws, f"{columns['purification']}{row}", None)
        if columns["scale"]:
            _set_cell_or_merged_anchor(ws, f"{columns['scale']}{row}", None)
        for key in ("mod5", "mod3", "dual_mod", "type", "remarks"):
            _set_cell_or_merged_anchor(ws, f"{columns[key]}{row}", None)

    for i, (pname, seq, warn) in enumerate(flat, start=1):
        row = START_ROW + (i - 1)
        _set_cell_or_merged_anchor(ws, f"{columns['id']}{row}", f"{i:02d}")
        _set_cell_or_merged_anchor(ws, f"{columns['name']}{row}", pname)
        _set_cell_or_merged_anchor(ws, f"{columns['sequence']}{row}", seq)
        _set_cell_or_merged_anchor(ws, f"{columns['length']}{row}", _len_formula(row))
        _set_cell_or_merged_anchor(ws, f"{columns['tubes']}{row}", tubes)
        if columns["scale"]:
            _set_cell_or_merged_anchor(ws, f"{columns['scale']}{row}", nmoles)
        else:
            _set_cell_or_merged_anchor(ws, f"{columns['amount']}{row}", nmoles)
        if columns["unit"]:
            _set_cell_or_merged_anchor(ws, f"{columns['unit']}{row}", "nmol")
        _set_cell_or_merged_anchor(ws, f"{columns['purification']}{row}", purification)
        _set_cell_or_merged_anchor(ws, f"{columns['remarks']}{row}", warn)

    _set_cell_or_merged_anchor(ws, "E16", f"=SUM(D{START_ROW}:D{needed_max_row})")

    bio = io.BytesIO()
    wb.save(bio)
    bio.seek(0)
    return bio.getvalue()


def export_sangon_template_bytes(
    primer_rows,
    purification: str = "HAP",
    nmoles: int = 5,
    tubes: int = 1,
    profile: dict | None = None,
    delivery_format: str | None = None,
    blank_profile_fields: bool = False,
):
    wb = openpyxl.load_workbook(TEMPLATE_PATH)
    return _export_with_workbook(
        wb,
        primer_rows=primer_rows,
        purification=purification,
        nmoles=nmoles,
        tubes=tubes,
        profile=profile,
        delivery_format=delivery_format,
        blank_profile_fields=blank_profile_fields,
    )


def export_from_template_bytes(
    primer_rows,
    template_bytes: bytes,
    purification: str = "HAP",
    nmoles: int = 5,
    tubes: int = 1,
    profile: dict | None = None,
    delivery_format: str | None = None,
    blank_profile_fields: bool = False,
):
    wb = openpyxl.load_workbook(io.BytesIO(template_bytes))
    return _export_with_workbook(
        wb,
        primer_rows=primer_rows,
        purification=purification,
        nmoles=nmoles,
        tubes=tubes,
        profile=profile,
        delivery_format=delivery_format,
        blank_profile_fields=blank_profile_fields,
    )
