import os
import sqlite3
from pathlib import Path
from config import load_config


def _resolve_db_path() -> Path:
    cfg = load_config()
    raw = (cfg.db_path or "").strip() or (os.environ.get("DB_PATH") or "").strip()
    if raw:
        p = Path(raw).expanduser()
        if p.suffix.lower() != ".db":
            p = p / "app.db"
        p.parent.mkdir(parents=True, exist_ok=True)
        return p
    return Path(__file__).with_name("app.db")


DB_PATH = _resolve_db_path()


def get_conn():
    conn = sqlite3.connect(DB_PATH, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys=ON")
    init_db(conn)
    conn.execute("PRAGMA foreign_keys=ON")
    return conn


def get_app_setting(conn: sqlite3.Connection, key: str, default: str | None = None) -> str | None:
    row = conn.execute("SELECT value FROM app_settings WHERE key=?", (key,)).fetchone()
    if not row:
        return default
    return row["value"]


def set_app_setting(conn: sqlite3.Connection, key: str, value: str | None):
    if value is None:
        conn.execute("DELETE FROM app_settings WHERE key=?", (key,))
    else:
        conn.execute(
            """
            INSERT INTO app_settings(key, value)
            VALUES(?, ?)
            ON CONFLICT(key) DO UPDATE SET value=excluded.value
            """,
            (key, value),
        )
    conn.commit()


def can_delete_gene_ids(conn: sqlite3.Connection, user_id: int, gene_ids: list[int]) -> tuple[bool, str, dict[int, list[int]]]:
    if not gene_ids:
        return True, "", {}
    placeholders = ",".join("?" for _ in gene_ids)
    sql = f"""
        SELECT g.id AS gene_id, p.run_id
        FROM genes g
        JOIN primers p ON p.gene_id = g.id
        JOIN primer_runs r ON r.id = p.run_id
        WHERE g.user_id=? AND r.user_id=? AND g.id IN ({placeholders})
        ORDER BY p.run_id DESC
    """
    rows = conn.execute(sql, [user_id, user_id, *gene_ids]).fetchall()
    if not rows:
        return True, "", {}

    refs: dict[int, list[int]] = {}
    for r in rows:
        gid = int(r["gene_id"])
        rid = int(r["run_id"])
        refs.setdefault(gid, [])
        if rid not in refs[gid]:
            refs[gid].append(rid)
    return False, "Selected genes are referenced by history runs.", refs


def can_delete_backbone(conn: sqlite3.Connection, user_id: int, backbone_id: int) -> tuple[bool, str, list[int]]:
    rows = conn.execute(
        """
        SELECT id
        FROM primer_runs
        WHERE user_id=? AND backbone_id=?
        ORDER BY id DESC
        """,
        (user_id, backbone_id),
    ).fetchall()
    if not rows:
        return True, "", []
    run_ids = [int(r["id"]) for r in rows]
    return False, "Backbone is referenced by history runs.", run_ids


def _backbones_schema_sql(conn: sqlite3.Connection) -> str:
    row = conn.execute("SELECT sql FROM sqlite_master WHERE type='table' AND name='backbones'").fetchone()
    return (row["sql"] or "") if row else ""


def _migrate_backbones_unique(conn: sqlite3.Connection):
    """
    Upgrade old backbones uniqueness:
      UNIQUE(user_id, name) -> UNIQUE(user_id, name, enzyme_sites)
    Keep backbone ids stable so primer_runs.backbone_id remains valid.
    """
    sql = _backbones_schema_sql(conn)
    if not sql:
        return
    normalized = " ".join(sql.replace("\n", " ").split()).lower()
    if "unique(user_id, name, enzyme_sites)" in normalized:
        return

    fk_before = int(conn.execute("PRAGMA foreign_keys").fetchone()[0])
    conn.execute("PRAGMA foreign_keys=OFF")
    cur = conn.cursor()
    cur.execute("BEGIN")
    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS backbones_new(
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            user_id INTEGER NOT NULL,
            name TEXT NOT NULL,
            enzyme_sites TEXT NOT NULL DEFAULT '',
            homo_f TEXT NOT NULL,
            homo_r TEXT NOT NULL,
            defaults_json TEXT,
            created_at TEXT NOT NULL,
            UNIQUE(user_id, name, enzyme_sites),
            FOREIGN KEY(user_id) REFERENCES users(id) ON DELETE CASCADE
        );
        """
    )
    cur.execute(
        """
        INSERT INTO backbones_new(id, user_id, name, enzyme_sites, homo_f, homo_r, defaults_json, created_at)
        SELECT id, user_id, name, COALESCE(enzyme_sites, ''), homo_f, homo_r, defaults_json, created_at
        FROM backbones;
        """
    )
    cur.execute("DROP TABLE backbones;")
    cur.execute("ALTER TABLE backbones_new RENAME TO backbones;")
    cur.execute("COMMIT")
    if fk_before:
        conn.execute("PRAGMA foreign_keys=ON")


def init_db(conn: sqlite3.Connection):
    cur = conn.cursor()

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS users(
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            username TEXT UNIQUE NOT NULL,
            password_hash TEXT NOT NULL,
            created_at TEXT NOT NULL
        );
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS user_profiles(
            user_id INTEGER PRIMARY KEY,
            orderer_name TEXT,
            mobile TEXT,
            email TEXT,
            legal_person TEXT,
            phone TEXT,
            institution TEXT,
            invoice_title TEXT,
            payment_method TEXT,
            invoice_with_goods TEXT,
            shipping_address TEXT,
            header_note TEXT,
            updated_at TEXT NOT NULL,
            FOREIGN KEY(user_id) REFERENCES users(id) ON DELETE CASCADE
        );
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS backbones(
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            user_id INTEGER NOT NULL,
            name TEXT NOT NULL,
            enzyme_sites TEXT NOT NULL DEFAULT '',
            homo_f TEXT NOT NULL,
            homo_r TEXT NOT NULL,
            defaults_json TEXT,
            created_at TEXT NOT NULL,
            UNIQUE(user_id, name, enzyme_sites),
            FOREIGN KEY(user_id) REFERENCES users(id) ON DELETE CASCADE
        );
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS genes(
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            user_id INTEGER NOT NULL,
            gene_name TEXT NOT NULL,
            sequence TEXT NOT NULL,
            created_at TEXT NOT NULL,
            UNIQUE(user_id, gene_name),
            FOREIGN KEY(user_id) REFERENCES users(id) ON DELETE CASCADE
        );
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS primer_runs(
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            user_id INTEGER NOT NULL,
            backbone_id INTEGER NOT NULL,
            params_json TEXT NOT NULL,
            created_at TEXT NOT NULL,
            note TEXT,
            FOREIGN KEY(user_id) REFERENCES users(id) ON DELETE CASCADE,
            FOREIGN KEY(backbone_id) REFERENCES backbones(id) ON DELETE CASCADE
        );
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS primers(
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            run_id INTEGER NOT NULL,
            gene_id INTEGER NOT NULL,
            primer_name_f TEXT NOT NULL,
            primer_name_r TEXT NOT NULL,
            gene_bind_f TEXT NOT NULL,
            gene_bind_r TEXT NOT NULL,
            final_f TEXT NOT NULL,
            final_r TEXT NOT NULL,
            bind_start_f INTEGER,
            bind_end_f INTEGER,
            bind_start_r INTEGER,
            bind_end_r INTEGER,
            amplicon_seq TEXT,
            tm_f REAL,
            tm_r REAL,
            gc_f REAL,
            gc_r REAL,
            amplicon_len INTEGER,
            warnings TEXT,
            selected_reason TEXT,
            tm_delta REAL,
            rescue_tier TEXT,
            tm_profile TEXT,
            tm_calc_params TEXT,
            FOREIGN KEY(run_id) REFERENCES primer_runs(id) ON DELETE CASCADE,
            FOREIGN KEY(gene_id) REFERENCES genes(id) ON DELETE CASCADE
        );
        """
    )

    primer_cols = {r[1] for r in cur.execute("PRAGMA table_info(primers)").fetchall()}
    if "amplicon_len" not in primer_cols:
        cur.execute("ALTER TABLE primers ADD COLUMN amplicon_len INTEGER")
    if "selected_reason" not in primer_cols:
        cur.execute("ALTER TABLE primers ADD COLUMN selected_reason TEXT")
    if "tm_delta" not in primer_cols:
        cur.execute("ALTER TABLE primers ADD COLUMN tm_delta REAL")
    if "rescue_tier" not in primer_cols:
        cur.execute("ALTER TABLE primers ADD COLUMN rescue_tier TEXT")
    if "tm_profile" not in primer_cols:
        cur.execute("ALTER TABLE primers ADD COLUMN tm_profile TEXT")
    if "tm_calc_params" not in primer_cols:
        cur.execute("ALTER TABLE primers ADD COLUMN tm_calc_params TEXT")
    if "bind_start_f" not in primer_cols:
        cur.execute("ALTER TABLE primers ADD COLUMN bind_start_f INTEGER")
    if "bind_end_f" not in primer_cols:
        cur.execute("ALTER TABLE primers ADD COLUMN bind_end_f INTEGER")
    if "bind_start_r" not in primer_cols:
        cur.execute("ALTER TABLE primers ADD COLUMN bind_start_r INTEGER")
    if "bind_end_r" not in primer_cols:
        cur.execute("ALTER TABLE primers ADD COLUMN bind_end_r INTEGER")
    if "amplicon_seq" not in primer_cols:
        cur.execute("ALTER TABLE primers ADD COLUMN amplicon_seq TEXT")

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS app_settings(
            key TEXT PRIMARY KEY,
            value TEXT
        );
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS order_templates(
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            user_id INTEGER NOT NULL,
            template_name TEXT NOT NULL,
            file_name TEXT NOT NULL,
            template_blob BLOB NOT NULL,
            is_default INTEGER NOT NULL DEFAULT 0,
            created_at TEXT NOT NULL,
            updated_at TEXT NOT NULL,
            UNIQUE(user_id, template_name),
            FOREIGN KEY(user_id) REFERENCES users(id) ON DELETE CASCADE
        );
        """
    )
    cur.execute(
        """
        CREATE INDEX IF NOT EXISTS idx_order_templates_user_default
        ON order_templates(user_id, is_default, id DESC);
        """
    )

    _migrate_backbones_unique(conn)
    conn.commit()
