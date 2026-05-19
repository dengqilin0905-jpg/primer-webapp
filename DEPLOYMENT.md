# Deployment Guide

## 1) Local app usage

### Prepare environment

```bash
python -m venv .venv
```

Windows:

```bash
.venv\Scripts\activate
```

Linux/macOS:

```bash
source .venv/bin/activate
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Run:

```bash
streamlit run app/app.py --server.address 127.0.0.1 --server.port 8501
```

## 2) GitHub update flow

Recommended workflow:

```bash
git checkout main
git pull origin main
git checkout -b feature/<your-feature>
# code + test
git add .
git commit -m "feat: <summary>"
git push -u origin feature/<your-feature>
```

Then open a PR into `main`.

## 3) GitHub Actions auto deploy to server

This repository supports deployment from GitHub Actions to your Linux server over SSH.

### Trigger

- Auto trigger: push to `main`
- Manual trigger: `workflow_dispatch`

Workflow file:

- `.github/workflows/deploy.yml`

### Required GitHub Secrets

Set these in **Repository Settings -> Secrets and variables -> Actions**:

- `SSH_HOST` (server public IP/domain)
- `SSH_USER` (SSH login user)
- `SSH_PORT` (usually `22`)
- `SSH_PRIVATE_KEY` (private key content, multi-line)

### Server-side behavior

GitHub Actions will:

1. Package current repo on runner (not on server).
2. Upload archive to `/tmp/primer-webapp-release/release.tar.gz`.
3. Upload and execute `scripts/remote_deploy.sh`.

`scripts/remote_deploy.sh` will:

1. stop service: `primer-webapp`
2. backup `/opt/primer-webapp` to `/opt/primer-webapp.bak_<timestamp>`
3. preserve `/opt/primer-webapp/.env.server`
4. preserve `/opt/primer-webapp/data/`
5. replace `/opt/primer-webapp` with uploaded release code
6. restore `.env.server` and `data/`
7. create/update `.venv` and install `requirements.txt`
8. set executable permissions on run scripts
9. `systemctl daemon-reload`
10. `systemctl start primer-webapp`
11. print `systemctl status primer-webapp --no-pager`

### Exclusion policy during packaging

The workflow excludes:

- `.git`
- `.venv`
- `__pycache__`
- `*.db`
- `.env*`
- `data/`
- `app/_deps/`
- `python/`
- temp/cache files

So existing server database and user data are preserved by restore logic.

## 4) Streamlit Community Cloud deployment

1. Push code to GitHub.
2. In Streamlit Community Cloud:
   - New app
   - Select repository and branch (recommend `main`)
   - Main file path: `app/app.py`
3. Add secrets/environment if needed (optional):
   - `APP_MODE`
   - `DB_PATH`
   - `APP_PASSWORD`
4. Deploy.

Notes:

- This project uses SQLite by default; for Streamlit Cloud, DB persistence depends on platform storage behavior.
- For lightweight internal demo it works; for larger multi-user production use, migrate DB backend.

## 5) Netlify statement

Netlify is not recommended for this current stack because:

- This app is Streamlit + SQLite (stateful Python server app).
- Netlify is better for static frontend or separated frontend/backend architecture.

## 6) Future public architecture recommendation

If you later build a public-scale website, prefer:

- Netlify frontend
- FastAPI backend
- PostgreSQL (Supabase / Neon / managed PG)

That split provides better scalability, persistence, and access control.

