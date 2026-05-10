import io
import time
import sqlite3
import streamlit as st
import openpyxl

from excel_export import TEMPLATE_PATH
from file_safety import validate_upload


ALLOWED_TEMPLATE_EXTS = {".xlsx"}


def _normalize_template_name(name: str) -> str:
    v = (name or "").strip()
    return v or "My template"


def _validate_template_bytes(data: bytes) -> None:
    try:
        wb = openpyxl.load_workbook(io.BytesIO(data))
    except Exception as exc:
        raise ValueError(f"Invalid Excel template: {exc}") from exc
    if not wb.sheetnames:
        raise ValueError("Template has no worksheet.")


def get_user_templates(conn: sqlite3.Connection, user_id: int):
    return conn.execute(
        """
        SELECT id, user_id, template_name, file_name, is_default, created_at, updated_at
        FROM order_templates
        WHERE user_id=?
        ORDER BY is_default DESC, updated_at DESC, id DESC
        """,
        (user_id,),
    ).fetchall()


def get_default_template_id(conn: sqlite3.Connection, user_id: int) -> int | None:
    row = conn.execute(
        "SELECT id FROM order_templates WHERE user_id=? AND is_default=1 ORDER BY id DESC LIMIT 1",
        (user_id,),
    ).fetchone()
    return int(row["id"]) if row else None


def get_template_bytes_by_id(conn: sqlite3.Connection, user_id: int, template_id: int) -> bytes | None:
    row = conn.execute(
        "SELECT template_blob FROM order_templates WHERE user_id=? AND id=?",
        (user_id, int(template_id)),
    ).fetchone()
    if not row:
        return None
    return bytes(row["template_blob"])


def get_effective_template_bytes(
    conn: sqlite3.Connection, user_id: int, template_id: int | None
) -> tuple[bytes, str]:
    if template_id:
        data = get_template_bytes_by_id(conn, user_id, int(template_id))
        if data:
            row = conn.execute(
                "SELECT template_name FROM order_templates WHERE user_id=? AND id=?",
                (user_id, int(template_id)),
            ).fetchone()
            name = row["template_name"] if row else f"Template#{template_id}"
            return data, name
    default_id = get_default_template_id(conn, user_id)
    if default_id:
        data = get_template_bytes_by_id(conn, user_id, default_id)
        if data:
            row = conn.execute(
                "SELECT template_name FROM order_templates WHERE user_id=? AND id=?",
                (user_id, default_id),
            ).fetchone()
            return data, (row["template_name"] if row else "Default template")
    return TEMPLATE_PATH.read_bytes(), "Built-in template"


def _set_default_template(conn: sqlite3.Connection, user_id: int, template_id: int) -> None:
    cur = conn.cursor()
    cur.execute("UPDATE order_templates SET is_default=0 WHERE user_id=?", (user_id,))
    cur.execute(
        "UPDATE order_templates SET is_default=1, updated_at=? WHERE user_id=? AND id=?",
        (time.strftime("%Y-%m-%d %H:%M:%S"), user_id, template_id),
    )
    conn.commit()


def templates_page(user, config):
    st.subheader("Order Templates")
    st.caption("Upload and manage your own Excel ordering templates.")

    from db import get_conn

    conn = get_conn()
    user_id = int(user["id"])

    quota = int(config.max_templates_per_user)
    size_limit = int(config.max_template_size_mb)

    with st.expander("Upload template", expanded=True):
        if quota > 0:
            st.caption(f"Per-user template quota: {quota}. File size limit: {size_limit} MB.")
        else:
            st.caption(f"No hard template count quota in local mode. File size limit: {size_limit} MB.")
        up = st.file_uploader(
            "Upload .xlsx template",
            type=["xlsx"],
            accept_multiple_files=False,
            key="template_upload_file",
        )
        tname = st.text_input("Template name", key="template_upload_name")
        make_default = st.checkbox("Set as default template", value=False, key="template_upload_default")
        if st.button("Save template", key="template_upload_save"):
            if not up:
                st.warning("Please choose an .xlsx file.")
            else:
                try:
                    validate_upload(up, ALLOWED_TEMPLATE_EXTS, size_limit)
                    blob = up.getvalue()
                    _validate_template_bytes(blob)
                    name = _normalize_template_name(tname or up.name.rsplit(".", 1)[0])
                    existing = conn.execute(
                        "SELECT id FROM order_templates WHERE user_id=? AND template_name=?",
                        (user_id, name),
                    ).fetchone()
                    count_row = conn.execute(
                        "SELECT COUNT(*) AS c FROM order_templates WHERE user_id=?",
                        (user_id,),
                    ).fetchone()
                    count = int(count_row["c"]) if count_row else 0
                    if quota > 0 and not existing and count >= quota:
                        raise ValueError(f"Template quota reached ({quota}). Delete one before uploading.")

                    now = time.strftime("%Y-%m-%d %H:%M:%S")
                    if existing:
                        conn.execute(
                            """
                            UPDATE order_templates
                            SET file_name=?, template_blob=?, updated_at=?
                            WHERE user_id=? AND id=?
                            """,
                            (up.name, blob, now, user_id, int(existing["id"])),
                        )
                        target_id = int(existing["id"])
                    else:
                        conn.execute(
                            """
                            INSERT INTO order_templates(user_id, template_name, file_name, template_blob, is_default, created_at, updated_at)
                            VALUES(?,?,?,?,0,?,?)
                            """,
                            (user_id, name, up.name, blob, now, now),
                        )
                        target_id = int(conn.execute("SELECT last_insert_rowid()").fetchone()[0])
                    conn.commit()

                    if make_default:
                        _set_default_template(conn, user_id, target_id)
                    st.success(f"Template saved: {name}")
                    st.rerun()
                except Exception as exc:
                    st.error(f"Failed to save template: {exc}")

    rows = get_user_templates(conn, user_id)
    if not rows:
        st.info("No user templates yet. Built-in Sangon template will be used.")
        conn.close()
        return

    st.markdown("#### Your templates")
    st.dataframe([dict(r) for r in rows], use_container_width=True)
    labels = [f"{r['id']} - {r['template_name']}" + (" (default)" if int(r["is_default"]) == 1 else "") for r in rows]
    selected_label = st.selectbox("Select template", labels, key="template_pick")
    selected_idx = labels.index(selected_label)
    selected = rows[selected_idx]
    selected_id = int(selected["id"])

    c1, c2, c3 = st.columns([1, 1, 1])
    with c1:
        if st.button("Set default", key="template_set_default"):
            _set_default_template(conn, user_id, selected_id)
            st.success("Default template updated.")
            st.rerun()
    with c2:
        data = get_template_bytes_by_id(conn, user_id, selected_id)
        if data:
            st.download_button(
                "Download template",
                data=data,
                file_name=selected["file_name"] or f"template_{selected_id}.xlsx",
                key=f"template_download_{selected_id}",
            )
    with c3:
        confirm_delete = st.checkbox("Confirm delete", key=f"template_confirm_delete_{selected_id}")
        if st.button("Delete template", key="template_delete_btn"):
            if not confirm_delete:
                st.warning("Please confirm deletion.")
            else:
                conn.execute("DELETE FROM order_templates WHERE user_id=? AND id=?", (user_id, selected_id))
                conn.commit()
                st.success("Template deleted.")
                st.rerun()
    conn.close()

