#!/usr/bin/env bash
set -euo pipefail

export APP_MODE="${APP_MODE:-server}"
export APP_DIR="${APP_DIR:-/opt/primer-webapp}"
export DB_PATH="${DB_PATH:-/opt/primer-webapp/data/app.db}"
export APP_HOST="${APP_HOST:-127.0.0.1}"
export APP_PORT="${APP_PORT:-8501}"
export MAX_UPLOAD_SIZE_MB="${MAX_UPLOAD_SIZE_MB:-20}"
export MAX_TEMPLATE_SIZE_MB="${MAX_TEMPLATE_SIZE_MB:-2}"
export MAX_TEMPLATES_PER_USER="${MAX_TEMPLATES_PER_USER:-5}"
export TEMP_FILE_RETENTION_HOURS="${TEMP_FILE_RETENTION_HOURS:-24}"
export TEMP_DIR="${TEMP_DIR:-/tmp/primer-webapp}"
export STREAMLIT_BIN="${STREAMLIT_BIN:-${APP_DIR}/.venv/bin/streamlit}"

mkdir -p "$(dirname "$DB_PATH")"
mkdir -p "$TEMP_DIR"

if [[ ! -x "$STREAMLIT_BIN" ]]; then
  echo "ERROR: streamlit executable not found: $STREAMLIT_BIN" >&2
  exit 1
fi

exec "$STREAMLIT_BIN" run app/app.py \
  --server.address "$APP_HOST" \
  --server.port "$APP_PORT" \
  --server.baseUrlPath primer \
  --server.maxUploadSize "$MAX_UPLOAD_SIZE_MB"
