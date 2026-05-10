# Primer WebApp

Streamlit-based primer design tool with user login, history tracking, multi-format sequence import, and Excel/TXT/FASTA export.

This repository now supports two operation modes:

- `local` mode: personal long-term local use.
- `server` mode: small-group cloud trial with stronger resource and isolation controls.

## Core Stack

- Language: Python
- Framework: Streamlit
- Database: SQLite
- Main entry: `app/app.py`

## Features

- User registration/login (password hash: PBKDF2-HMAC-SHA256 + random salt)
- Backbone management (supports same backbone name with different enzyme sites)
- Gene import: FASTA/TXT/GenBank/SnapGene `.dna` (with content sniffing)
- Mapping-based batch primer design (`CDS` / `Promoter` / `Other`)
- Per-mapping `remove_stop` for CDS
- History with re-export
- Export:
  - Excel order sheet
  - Fragment TXT (full/simple)
  - non-CDS amplicon FASTA
- User template system:
  - upload multiple Excel templates
  - set default template
  - delete/download template
  - choose template when exporting (design + history)

## Repository Layout

```text
app/
  app.py
  db.py
  primer_design.py
  design_ui.py
  history_ui.py
  templates_ui.py
  ...
scripts/
  run_local.sh
  run_server.sh
deploy/
  systemd/primer-webapp.service
  nginx/primer-webapp.conf
examples/
docs/
```

## Local Mode (GitHub download, run on your own machine)

### 1) Create venv and install

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### 2) Start app

```bash
export APP_MODE=local
export DB_PATH=./data/app.db
streamlit run app/app.py --server.address 127.0.0.1 --server.port 8501
```

Open `http://127.0.0.1:8501`.

## Server Mode (Alibaba Cloud Ubuntu 22.04)

### 1) Server requirements

- Ubuntu 22.04
- Python 3.10+
- nginx + systemd

### 2) Deploy path example

```bash
sudo mkdir -p /opt/primer-webapp
sudo chown -R $USER:$USER /opt/primer-webapp
cd /opt/primer-webapp
git clone <your-repo-url> .
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
chmod +x scripts/run_local.sh scripts/run_server.sh
```

### 3) Environment (server mode)

```bash
export APP_MODE=server
export DB_PATH=/opt/primer-webapp/data/app.db
export APP_HOST=127.0.0.1
export APP_PORT=8501
export MAX_TEMPLATE_SIZE_MB=2
export MAX_TEMPLATES_PER_USER=5
export MAX_UPLOAD_SIZE_MB=20
export TEMP_FILE_RETENTION_HOURS=24
export TEMP_DIR=/tmp/primer-webapp
# optional shared gate
export APP_PASSWORD=<optional-site-password>
```

### 4) Start with systemd

Use `deploy/systemd/primer-webapp.service` and set:

- `User`
- `WorkingDirectory`
- `EnvironmentFile` path

Then:

```bash
sudo systemctl daemon-reload
sudo systemctl enable primer-webapp
sudo systemctl start primer-webapp
sudo systemctl status primer-webapp
```

### 5) nginx reverse proxy

Use `deploy/nginx/primer-webapp.conf`, then:

```bash
sudo ln -s /etc/nginx/sites-available/primer-webapp.conf /etc/nginx/sites-enabled/primer-webapp.conf
sudo nginx -t
sudo systemctl reload nginx
```

## Security and Isolation Notes

- Passwords are not stored in plaintext.
- DB operations for history/templates/genes/backbones are user-scoped.
- Server upload checks:
  - extension allowlist
  - file size limit (`MAX_UPLOAD_SIZE_MB`, `MAX_TEMPLATE_SIZE_MB`)
- Server template quota per user: `MAX_TEMPLATES_PER_USER` (default 5)
- Template data is bound to `user_id`.
- Temporary files are isolated in `TEMP_DIR` and cleaned by retention policy.
- Keep primer app and private JBrowse in separate directories and separate nginx site configs.

## Config Variables

- `APP_MODE` = `local` or `server`
- `DB_PATH`
- `APP_HOST`
- `APP_PORT`
- `MAX_TEMPLATE_SIZE_MB`
- `MAX_TEMPLATES_PER_USER`
- `MAX_UPLOAD_SIZE_MB`
- `TEMP_FILE_RETENTION_HOURS`
- `TEMP_DIR`
- `APP_PASSWORD` (optional)

## Input / Output

- Inputs:
  - gene files (`.fasta/.fa/.fna/.txt/.gb/.gbk/.genbank/.dna`)
  - backbone settings
  - design parameters
- Outputs:
  - Excel order sheet (`.xlsx`, selected template)
  - fragment record TXT
  - non-CDS amplicon FASTA

## Release Suggestion

Recommend publishing `v1.0.0` with:

1. source code
2. local setup guide
3. server deployment guide
4. examples
5. changelog
