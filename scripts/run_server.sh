#!/usr/bin/env bash
set -euo pipefail

export APP_MODE="${APP_MODE:-server}"
export DB_PATH="${DB_PATH:-/opt/primer-webapp/data/app.db}"
export APP_HOST="${APP_HOST:-127.0.0.1}"
export APP_PORT="${APP_PORT:-8501}"
export MAX_UPLOAD_SIZE_MB="${MAX_UPLOAD_SIZE_MB:-20}"
export MAX_TEMPLATE_SIZE_MB="${MAX_TEMPLATE_SIZE_MB:-2}"
export MAX_TEMPLATES_PER_USER="${MAX_TEMPLATES_PER_USER:-5}"
export TEMP_FILE_RETENTION_HOURS="${TEMP_FILE_RETENTION_HOURS:-24}"
export TEMP_DIR="${TEMP_DIR:-/tmp/primer-webapp}"

mkdir -p "$(dirname "$DB_PATH")"
mkdir -p "$TEMP_DIR"

exec streamlit run app/app.py \
  --server.address "$APP_HOST" \
  --server.port "$APP_PORT" \
  --server.maxUploadSize "$MAX_UPLOAD_SIZE_MB"
