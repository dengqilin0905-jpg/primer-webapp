import os
import sqlite3
from pathlib import Path

def _resolve_db_path() -> Path:
    # Prefer env var (for deployment). Accept either a file path or a directory.
    raw = (os.environ.get("DB_PATH") or "").strip()
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
    init_db(conn)
    return conn

def init_db(conn: sqlite3.Connection):
    cur = conn.cursor()

    cur.execute('''
    CREATE TABLE IF NOT EXISTS users(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        username TEXT UNIQUE NOT NULL,
        password_hash TEXT NOT NULL,
        created_at TEXT NOT NULL
    );
    ''')

    cur.execute('''
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
    ''')

    cur.execute('''
    CREATE TABLE IF NOT EXISTS backbones(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        user_id INTEGER NOT NULL,
        name TEXT NOT NULL,
        enzyme_sites TEXT,
        homo_f TEXT NOT NULL,
        homo_r TEXT NOT NULL,
        defaults_json TEXT,
        created_at TEXT NOT NULL,
        UNIQUE(user_id, name),
        FOREIGN KEY(user_id) REFERENCES users(id) ON DELETE CASCADE
    );
    ''')

    cur.execute('''
    CREATE TABLE IF NOT EXISTS genes(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        user_id INTEGER NOT NULL,
        gene_name TEXT NOT NULL,
        sequence TEXT NOT NULL,
        created_at TEXT NOT NULL,
        UNIQUE(user_id, gene_name),
        FOREIGN KEY(user_id) REFERENCES users(id) ON DELETE CASCADE
    );
    ''')

    cur.execute('''
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
    ''')

    cur.execute('''
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
        tm_f REAL,
        tm_r REAL,
        gc_f REAL,
        gc_r REAL,
        warnings TEXT,
        FOREIGN KEY(run_id) REFERENCES primer_runs(id) ON DELETE CASCADE,
        FOREIGN KEY(gene_id) REFERENCES genes(id) ON DELETE CASCADE
    );
    ''')

    conn.commit()
