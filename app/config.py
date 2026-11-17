import os
from dataclasses import dataclass


def _env_str(key: str, default: str) -> str:
    val = (os.environ.get(key) or "").strip()
    return val if val else default


def _env_int(key: str, default: int, min_value: int | None = None) -> int:
    raw = (os.environ.get(key) or "").strip()
    if not raw:
        return default
    try:
        value = int(raw)
    except ValueError:
        return default
    if min_value is not None and value < min_value:
        return min_value
    return value


@dataclass(frozen=True)
class AppConfig:
    app_mode: str
    db_path: str
    app_host: str
    app_port: int
    max_template_size_mb: int
    max_templates_per_user: int
    max_upload_size_mb: int
    temp_file_retention_hours: int
    temp_dir: str
    app_password: str

    @property
    def is_server(self) -> bool:
        return self.app_mode == "server"


def load_config() -> AppConfig:
    mode = _env_str("APP_MODE", "local").lower()
    if mode not in {"local", "server"}:
        mode = "local"

    if mode == "server":
        default_temp_dir = "/tmp/primer-webapp"
        default_template_quota = 5
        default_upload_mb = 20
    else:
        from resource_config import state_root
        default_temp_dir = str(state_root()/"tmp")
        default_template_quota = 0  # 0 means unlimited
        default_upload_mb = 50

    return AppConfig(
        app_mode=mode,
        db_path=_env_str("DB_PATH", ""),
        app_host=_env_str("APP_HOST", "127.0.0.1"),
        app_port=_env_int("APP_PORT", 8501, min_value=1),
        max_template_size_mb=_env_int("MAX_TEMPLATE_SIZE_MB", 2, min_value=1),
        max_templates_per_user=_env_int("MAX_TEMPLATES_PER_USER", default_template_quota, min_value=0),
        max_upload_size_mb=_env_int("MAX_UPLOAD_SIZE_MB", default_upload_mb, min_value=1),
        temp_file_retention_hours=_env_int("TEMP_FILE_RETENTION_HOURS", 24, min_value=1),
        temp_dir=_env_str("TEMP_DIR", default_temp_dir),
        app_password=_env_str("APP_PASSWORD", ""),
    )

