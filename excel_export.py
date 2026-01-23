import io
import re
from pathlib import Path
import openpyxl
from copy import copy

TEMPLATE_PATH = Path(__file__).with_name("sangon_template.xlsx")
SHEET_NAME = "引物合成订购表"
START_ROW = 18
BASE_TEMPLATE_MAX_ROW = 223  # original pre-formatted range

def _len_formula(row: int) -> str:
    return f"=LEN(C{row})-SUM((LEN(C{row})-LEN(SUBSTITUTE(C{row},修饰名称!E:E,修饰名称!F:F))))"

def _fill_header(ws, profile: dict | None):
    if not profile:
        return
    mapping = {
        'orderer_name': 'B4',
        'mobile': 'B5',
        'email': 'B6',
        'legal_person': 'B7',
        'phone': 'B8',
        'institution': 'B9',
        'invoice_title': 'B10',
        'payment_method': 'B11',
        'invoice_with_goods': 'B12',
        'shipping_address': 'B13',
        'header_note': 'B14',
    }
    for k, cell in mapping.items():
        v = profile.get(k)
        if v is not None and v != '':
            ws[cell].value = v

def _copy_row_style(ws, src_row: int, dst_row: int, max_col: int = 13):
    if ws.row_dimensions.get(src_row) and ws.row_dimensions[src_row].height:
        ws.row_dimensions[dst_row].height = ws.row_dimensions[src_row].height
    for c in range(1, max_col+1):
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
    for row in range(BASE_TEMPLATE_MAX_ROW+1, needed_max_row+1):
        _copy_row_style(ws, src, row)
        ws[f"D{row}"].value = _len_formula(row)

def export_sangon_template_bytes(primer_rows, purification: str = "HAP", nmoles: int = 5, tubes: int = 1, profile: dict | None = None):
    wb = openpyxl.load_workbook(TEMPLATE_PATH)
    ws = wb[SHEET_NAME]

    _fill_header(ws, profile)

    # Flatten: each primer row becomes two order rows (F then R)
    flat = []
    for r in primer_rows:
        warn = r.get("warnings","")
        flat.append((r["primer_name_f"], r["final_f"], warn))
        flat.append((r["primer_name_r"], r["final_r"], warn))

    needed_max_row = START_ROW + len(flat) - 1 if flat else START_ROW
    _ensure_rows(ws, needed_max_row)

    # Clear and fill
    for row in range(START_ROW, needed_max_row + 1):
        ws[f"A{row}"].value = None
        ws[f"B{row}"].value = None
        ws[f"C{row}"].value = None
        ws[f"D{row}"].value = _len_formula(row)
        ws[f"E{row}"].value = None
        ws[f"F{row}"].value = None
        ws[f"G{row}"].value = None
        ws[f"H{row}"].value = None
        ws[f"M{row}"].value = None

    for i, (pname, seq, warn) in enumerate(flat, start=1):
        row = START_ROW + (i - 1)
        ws[f"A{row}"].value = f"{i:02d}"
        ws[f"B{row}"].value = pname
        ws[f"C{row}"].value = seq
        ws[f"D{row}"].value = _len_formula(row)
        ws[f"E{row}"].value = tubes
        ws[f"G{row}"].value = purification
        ws[f"H{row}"].value = nmoles
        ws[f"M{row}"].value = warn

    try:
        ws["E16"].value = f"=SUM(D{START_ROW}:D{needed_max_row})"
    except Exception:
        pass

    bio = io.BytesIO()
    wb.save(bio)
    bio.seek(0)
    return bio.getvalue()
