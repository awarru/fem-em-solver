#!/usr/bin/env bash
set -euo pipefail

# One-command bootstrap + pending-test runner for fresh clones.
# Usage:
#   ./bootstrap_pending_tests.sh
#   ./bootstrap_pending_tests.sh --list
#   ./bootstrap_pending_tests.sh --chunk E3
#   ./bootstrap_pending_tests.sh --dry-run

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT_DIR"

COMPOSE_FILE_DEFAULT="$ROOT_DIR/docker/docker-compose.yml"
export COMPOSE_FILE="${COMPOSE_FILE:-$COMPOSE_FILE_DEFAULT}"

echo "[bootstrap] Repo root: $ROOT_DIR"
echo "[bootstrap] Compose file: $COMPOSE_FILE"

if ! command -v docker >/dev/null 2>&1; then
  echo "[bootstrap] ERROR: docker not found on PATH" >&2
  exit 2
fi

echo "[bootstrap] Ensuring fem-em-solver container is ready..."
docker compose up -d --build fem-em-solver

echo "[bootstrap] Running pending test flow..."
exec "$ROOT_DIR/run_tests.sh" "$@"
