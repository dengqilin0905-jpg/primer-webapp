# Changelog

## v1.1.1 - 2026-05-19

### Added

- GitHub Actions auto deploy workflow: `.github/workflows/deploy.yml`
  - Trigger on `push` to `main`
  - Uses SSH Secrets: `SSH_HOST`, `SSH_USER`, `SSH_PORT`, `SSH_PRIVATE_KEY`
  - Packages code on runner and deploys via SSH upload + remote script execution
- Remote deploy script: `scripts/remote_deploy.sh`
  - Stops `primer-webapp`
  - Backs up `/opt/primer-webapp`
  - Preserves/restores `.env.server` and `data/`
  - Rebuilds `.venv`, installs dependencies, restarts service

### Changed

- `DEPLOYMENT.md` expanded with GitHub Actions secrets setup and auto-deploy usage.

## v1.1.0 - 2026-05-18

### Added

- Design dual batch modes:
  - `mapping` (existing behavior)
  - `cartesian` (Gene × Backbone full combinations)
- Run status tracking in `primer_runs`: `pending/success/fail`.
- Backbone map fields:
  - `vector_sequence`
  - `insertion_mode`
  - `restriction_site_f/r`
  - `junction_extra_bases_f/r`
  - `vector_map_note`
- GenBank map export (`.gbk`) and multi-record ZIP export:
  - available in Design result and History re-export
  - SnapGene-compatible GenBank output with required features
- Optional frame-check settings in Design:
  - `enable_frame_check`
  - `construct_type`
  - warning integration into result/history/excel remarks
- Indentation checker script: `scripts/check_indent.py`.
- New docs: `DEPLOYMENT.md`.

### Changed

- Primer row schema extended to store per-result backbone/context data for mixed-backbone runs.
- Design run execution now creates one run per batch and auto-updates run status.
- History view now reads per-primer backbone info and supports GenBank map re-export.
- README rewritten with local usage, `.dna` import, batch design, and export guidance.

### Fixed

- Robust Tm formatting in History for `None/NaN` values.
- `.dna` parser path remains compatible with `snapgene-reader==0.1.23` and graceful file-level failure handling.
- User data isolation maintained in updated queries and exports.

## v1.0.0 - 2026-05-09

### Added

- Dual-mode configuration foundation (`local` / `server`) via `app/config.py`.
- File safety utilities (`app/file_safety.py`):
  - upload extension and size checks
  - temp directory management and cleanup
- User template system:
  - new DB table `order_templates`
  - template management UI (`app/templates_ui.py`)
  - per-user upload/default/delete/download
  - per-user template quota and size constraints (config-driven)
- Export supports custom template bytes:
  - `export_from_template_bytes(...)` in `app/excel_export.py`
- Template selection added to:
  - design-page export
  - history re-export
- Repository engineering files:
  - `README.md`
  - root `requirements.txt`
  - `.gitignore`
  - deployment scripts/config templates in `scripts/` and `deploy/`
  - `examples/` and `docs/`
  - `LICENSE`

### Changed

- `app/app.py` now:
  - loads central config
  - adds Templates tab
  - runs temp cleanup on startup in server mode
- `app/genes_ui.py` now enforces upload size/type checks.
- `app/history_ui.py` query tightened with explicit `user_id` scoping on primer rows.

### Compatibility

- SQLite migration is backward compatible:
  - existing databases auto-create new `order_templates` table if missing
  - no destructive schema rewrite for existing primer/backbone/history data

