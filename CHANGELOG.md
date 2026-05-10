# Changelog

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

