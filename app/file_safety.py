import os
import shutil
import tempfile
import time
import uuid
from pathlib import Path


def validate_upload(uploaded_file, allowed_ext: set[str], max_size_mb: int) -> None:
    name = (uploaded_file.name or "").strip()
    if "." not in name:
        raise ValueError("Uploaded file has no extension.")
    ext = "." + name.rsplit(".", 1)[-1].lower()
    if allowed_ext and ext not in allowed_ext:
        allowed = ", ".join(sorted(allowed_ext))
        raise ValueError(f"Unsupported file extension: {ext}. Allowed: {allowed}")
    size = getattr(uploaded_file, "size", None)
    if size is None:
        return
    max_bytes = int(max_size_mb * 1024 * 1024)
    if size > max_bytes:
        raise ValueError(f"File too large: {size} bytes > {max_bytes} bytes ({max_size_mb} MB).")


def ensure_temp_dir(path_hint: str) -> Path:
    raw = (path_hint or "").strip()
    p = Path(raw) if raw else Path(tempfile.gettempdir()) / "primer-webapp"
    p.mkdir(parents=True, exist_ok=True)
    return p


def make_temp_file_path(temp_dir: Path, suffix: str) -> Path:
    safe_suffix = suffix if suffix.startswith(".") else f".{suffix}"
    return temp_dir / f"{uuid.uuid4().hex}{safe_suffix}"


def cleanup_old_temp_files(temp_dir: Path, retention_hours: int) -> dict:
    now = time.time()
    ttl = max(1, int(retention_hours)) * 3600
    removed_files = 0
    removed_dirs = 0
    errors = 0
    if not temp_dir.exists():
        return {"removed_files": 0, "removed_dirs": 0, "errors": 0}

    for child in temp_dir.iterdir():
        try:
            mtime = child.stat().st_mtime
            if now - mtime <= ttl:
                continue
            if child.is_file():
                child.unlink(missing_ok=True)
                removed_files += 1
            elif child.is_dir():
                shutil.rmtree(child, ignore_errors=False)
                removed_dirs += 1
        except Exception:
            errors += 1
    return {"removed_files": removed_files, "removed_dirs": removed_dirs, "errors": errors}

