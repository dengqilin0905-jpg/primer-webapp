#!/usr/bin/env bash
set -euo pipefail

SERVICE_NAME="${SERVICE_NAME:-primer-webapp}"
APP_DIR="${APP_DIR:-/opt/primer-webapp}"
TMP_RELEASE_DIR="${TMP_RELEASE_DIR:-/tmp/primer-webapp-release}"
ARCHIVE_PATH="${1:-$TMP_RELEASE_DIR/release.tar.gz}"

if [[ ! -f "$ARCHIVE_PATH" ]]; then
  echo "ERROR: release archive not found: $ARCHIVE_PATH" >&2
  exit 1
fi

timestamp="$(date +%Y%m%d_%H%M%S)"
backup_dir="${APP_DIR}.bak_${timestamp}"
extract_dir="${TMP_RELEASE_DIR}/extract_${timestamp}"
restore_env="${TMP_RELEASE_DIR}/restore_env_${timestamp}"
restore_data="${TMP_RELEASE_DIR}/restore_data_${timestamp}"

mkdir -p "$TMP_RELEASE_DIR"
rm -rf "$extract_dir"
sudo rm -rf "$restore_env" "$restore_data"
mkdir -p "$extract_dir"

echo "[deploy] stop service: ${SERVICE_NAME}"
sudo systemctl stop "$SERVICE_NAME" || true

if [[ -d "$APP_DIR" ]]; then
  echo "[deploy] backup current app dir -> $backup_dir"
  sudo cp -a "$APP_DIR" "$backup_dir"
fi

if [[ -f "${APP_DIR}/.env.server" ]]; then
  echo "[deploy] preserve .env.server"
  sudo cp "${APP_DIR}/.env.server" "$restore_env"
fi

if [[ -d "${APP_DIR}/data" ]]; then
  echo "[deploy] preserve data directory"
  sudo cp -a "${APP_DIR}/data" "$restore_data"
fi

echo "[deploy] extract new release"
tar -xzf "$ARCHIVE_PATH" -C "$extract_dir"

echo "[deploy] replace app dir content"
sudo mkdir -p "$APP_DIR"
sudo find "$APP_DIR" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
sudo cp -a "${extract_dir}/." "$APP_DIR/"

if [[ -f "$restore_env" ]]; then
  echo "[deploy] restore .env.server"
  sudo cp "$restore_env" "${APP_DIR}/.env.server"
fi

if [[ -d "$restore_data" ]]; then
  echo "[deploy] restore data directory"
  sudo rm -rf "${APP_DIR}/data"
  sudo cp -a "$restore_data" "${APP_DIR}/data"
fi

echo "[deploy] create/update venv"
if [[ ! -d "${APP_DIR}/.venv" ]]; then
  python3 -m venv "${APP_DIR}/.venv"
fi
source "${APP_DIR}/.venv/bin/activate"
pip install --upgrade pip
pip install -r "${APP_DIR}/requirements.txt"

echo "[deploy] fix script executable bits"
chmod +x "${APP_DIR}/scripts/run_local.sh" "${APP_DIR}/scripts/run_server.sh" "${APP_DIR}/scripts/remote_deploy.sh"

echo "[deploy] restart service"
sudo systemctl daemon-reload
sudo systemctl start "$SERVICE_NAME"
sudo systemctl status "$SERVICE_NAME" --no-pager

echo "[deploy] done"
