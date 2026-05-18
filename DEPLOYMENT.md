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

## 3) Streamlit Community Cloud deployment

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

## 4) Netlify statement

Netlify is not recommended for this current stack because:

- This app is Streamlit + SQLite (stateful Python server app).
- Netlify is better for static frontend or separated frontend/backend architecture.

## 5) Future public architecture recommendation

If you later build a public-scale website, prefer:

- Netlify frontend
- FastAPI backend
- PostgreSQL (Supabase / Neon / managed PG)

That split provides better scalability, persistence, and access control.

