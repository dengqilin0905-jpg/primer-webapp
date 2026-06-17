# Changelog

## Unreleased

### Added

- Annotated construct map export now preserves readable backbone features from GenBank/SnapGene uploads when available.
- Optional SnapGene Desktop CLI conversion can export `.dna` maps from generated `.gbk` files.
- Added `scripts/qa_map_export.py` to validate primer generation, Sangon Excel export, GenBank map annotations, and optional SnapGene conversion with the canonical `test/` data.

### Changed

- Design and History map export controls now include insert feature type/name and advanced homology-arm annotation options.
- README, DEPLOYMENT, and AGENT rules now document local SnapGene Desktop requirements and the GitHub-vs-website deployment boundary.

### Fixed

- Sangon Excel exports now reset order header cells `B4:B14` to `0` before applying the current user's saved Profile fields, preventing personal header details from leaking through built-in or uploaded templates.

## v1.1.2 - 2026-05-26

### Fixed

- Sangon Excel export now clears the full order data area before writing current primers, preventing old template/order rows from leaking into new exports.
- CDS primer design now recovers from strict-filter failure on high-GC 5' starts by selecting the closest end-anchored pair with explicit warnings instead of failing the whole run.
- Primer design failure messages now include actionable Chinese guidance for adjusting Tm, primer length, GC range, homology arms, and stop-codon handling.

### Changed

- Backbones page now has Basic / Advanced setup modes.
- Basic backbone setup can auto-generate HOMO_F/HOMO_R from a vector sequence when the selected restriction sites are unique.
- Design page now keeps experimental decisions visible in Basic parameters and folds chemistry/frame-check controls into Advanced parameters.
- `AGENT.md` now documents required reproduction, test-folder, error-message, UI, export, and release rules.

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

