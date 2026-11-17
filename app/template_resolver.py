"""Shared order-template resolution for Web, batch, and external callers."""

from __future__ import annotations

import io
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from resource_config import resource_path,REPO_ROOT

import openpyxl


GLOBAL_DEFAULT_TEMPLATE_PATH = resource_path('order_template',REPO_ROOT/'resources/templates/order_generic.xlsx')


@dataclass(frozen=True)
class ResolvedOrderTemplate:
    template_bytes: bytes
    source_type: str
    template_name: str
    template_id: int | None = None
    resolved_path: Path | None = None


def _validate_xlsx_bytes(data: bytes, source: str) -> bytes:
    if not data:
        raise ValueError(f"Order template is empty: {source}")
    try:
        wb = openpyxl.load_workbook(io.BytesIO(data), read_only=False, data_only=False)
    except Exception as exc:
        raise ValueError(f"Invalid .xlsx order template ({source}): {exc}") from exc
    if not wb.sheetnames:
        raise ValueError(f"Order template has no worksheet: {source}")
    ws = wb[wb.sheetnames[0]]
    if ws.max_row < 18 or ws.max_column < 3:
        raise ValueError(f"Order template lacks the required order area: {source}")
    return data


def _read_template_path(path: Path) -> bytes:
    if not path.exists():
        raise FileNotFoundError(f"Order template path does not exist: {path}")
    if not path.is_file():
        raise ValueError(f"Order template path is not a file: {path}")
    if path.suffix.lower() != ".xlsx":
        raise ValueError(f"Order template must be an .xlsx file: {path}")
    return _validate_xlsx_bytes(path.read_bytes(), str(path))


def get_user_templates(conn: sqlite3.Connection, user_id: int):
    return conn.execute(
        """
        SELECT id, user_id, template_name, file_name, is_default, created_at, updated_at
        FROM order_templates
        WHERE user_id=?
        ORDER BY is_default DESC, updated_at DESC, id DESC
        """,
        (int(user_id),),
    ).fetchall()


def get_default_template_id(conn: sqlite3.Connection, user_id: int) -> int | None:
    row = conn.execute(
        "SELECT id FROM order_templates WHERE user_id=? AND is_default=1 ORDER BY id DESC LIMIT 1",
        (int(user_id),),
    ).fetchone()
    return int(row["id"]) if row else None


def _get_owned_template(conn: sqlite3.Connection, user_id: int, template_id: int):
    return conn.execute(
        """
        SELECT id, template_name, file_name, template_blob
        FROM order_templates
        WHERE user_id=? AND id=?
        """,
        (int(user_id), int(template_id)),
    ).fetchone()


def get_template_bytes_by_id(
    conn: sqlite3.Connection, user_id: int, template_id: int
) -> bytes | None:
    row = _get_owned_template(conn, user_id, template_id)
    return bytes(row["template_blob"]) if row else None


def resolve_order_template(
    *,
    template_path: str | Path | None = None,
    template_id: int | None = None,
    user_id: int | None = None,
    conn: sqlite3.Connection | None = None,
    db_path: str | Path | None = None,
) -> ResolvedOrderTemplate:
    """Resolve an order template without Streamlit session state."""
    owned_conn = None
    if db_path is not None and conn is None:
        owned_conn = sqlite3.connect(Path(db_path).expanduser())
        owned_conn.row_factory = sqlite3.Row
        conn = owned_conn

    try:
        return _resolve_order_template(
            template_path=template_path,
            template_id=template_id,
            user_id=user_id,
            conn=conn,
        )
    finally:
        if owned_conn is not None:
            owned_conn.close()


def _resolve_order_template(
    *,
    template_path: str | Path | None,
    template_id: int | None,
    user_id: int | None,
    conn: sqlite3.Connection | None,
) -> ResolvedOrderTemplate:
    if template_path is not None:
        path = Path(template_path).expanduser()
        data = _read_template_path(path)
        source = "global_default" if path == GLOBAL_DEFAULT_TEMPLATE_PATH else "explicit_path"
        return ResolvedOrderTemplate(data, source, path.stem, resolved_path=path)

    if template_id is not None:
        if user_id is None:
            raise ValueError("template_id requires user_id for ownership isolation")
        if conn is None:
            raise ValueError("template_id requires a database connection")
        row = _get_owned_template(conn, int(user_id), int(template_id))
        if not row:
            raise LookupError(
                f"Template {int(template_id)} does not belong to user {int(user_id)}"
            )
        data = _validate_xlsx_bytes(bytes(row["template_blob"]), f"template_id={int(template_id)}")
        return ResolvedOrderTemplate(
            data,
            "user_selected",
            row["template_name"] or row["file_name"] or f"Template#{int(template_id)}",
            template_id=int(row["id"]),
        )

    if user_id is not None and conn is not None:
        default_id = get_default_template_id(conn, int(user_id))
        if default_id is not None:
            row = _get_owned_template(conn, int(user_id), default_id)
            if row:
                data = _validate_xlsx_bytes(bytes(row["template_blob"]), f"template_id={default_id}")
                return ResolvedOrderTemplate(
                    data,
                    "user_default",
                    row["template_name"] or row["file_name"] or f"Template#{default_id}",
                    template_id=default_id,
                )

    data = _read_template_path(GLOBAL_DEFAULT_TEMPLATE_PATH)
    return ResolvedOrderTemplate(
        data,
        "global_default",
        GLOBAL_DEFAULT_TEMPLATE_PATH.stem,
        resolved_path=GLOBAL_DEFAULT_TEMPLATE_PATH,
    )
