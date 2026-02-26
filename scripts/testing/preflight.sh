#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
PENDING_FILE="$ROOT_DIR/docs/testing/pending-tests.md"
LOG_DIR="$ROOT_DIR/docs/testing/logs"
INDEX_FILE="$ROOT_DIR/docs/testing/test-results.md"
COMPOSE_FILE_DEFAULT="$ROOT_DIR/docker/docker-compose.yml"

mkdir -p "$LOG_DIR"

if [[ ! -f "$PENDING_FILE" ]]; then
  echo "Missing: $PENDING_FILE" >&2
  exit 2
fi

if [[ ! -f "$COMPOSE_FILE_DEFAULT" ]]; then
  echo "Missing docker compose file: $COMPOSE_FILE_DEFAULT" >&2
  exit 2
fi

if ! command -v docker >/dev/null 2>&1; then
  echo "docker is not installed or not on PATH" >&2
  exit 2
fi

if ! docker compose version >/dev/null 2>&1; then
  echo "docker compose is not available (docker compose version failed)" >&2
  exit 2
fi

if [[ ! -f "$INDEX_FILE" ]]; then
  cat > "$INDEX_FILE" <<'EOF'
# FEM-EM Manual Test Results

This file is append-only and records test runs executed by the human operator.
Each run has a full log in `docs/testing/logs/`.

| UTC Timestamp | Chunk | Exit | Log |
|---|---|---:|---|
EOF
fi

echo "Preflight OK"
echo "- Root: $ROOT_DIR"
echo "- Pending tests: $PENDING_FILE"
echo "- Logs dir: $LOG_DIR"
echo "- Index: $INDEX_FILE"
echo "- Compose file: $COMPOSE_FILE_DEFAULT"