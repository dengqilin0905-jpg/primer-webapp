# Primer WebApp (Streamlit + SQLite)

Features
- User login (local SQLite)
- Backbone library per user (homology arms)
- Batch gene import (FASTA/GenBank; .dna best-effort)
- Mapping-based combinations (Backbone → selected genes) — no forced Cartesian product
- Primer naming enforced: BackboneName_GeneName-F/R
- Export Sangon Excel using embedded template, with auto-extended rows beyond original formatting range
- Per-user Profile auto-fills template header cells B4–B14

Run
```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m streamlit run app.py
```


## Gene naming
Imported genes will be named as short aliases like OsTF1 (tail IDs are ignored). Use Utilities → Normalize existing gene names if needed.


## Deploy to Render (A2-1)
1) Push this folder to a GitHub repo.
2) In Render: New + → Web Service → connect repo.
3) Environment: Docker (auto-detected). Render will build using Dockerfile.
4) Add a Persistent Disk mounted at /data (1GB is enough) and set env var DB_PATH=/data/app.db.
5) Set APP_PASSWORD in Render (Environment → Secret) to restrict access.
