#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   scripts/testing/run_and_log.sh <chunk_id> <command...>
#   scripts/testing/run_and_log.sh --dry-run <chunk_id> <command...>
#   FEM_SOLVER_DRY_RUN=1 scripts/testing/run_and_log.sh <chunk_id> <command...>
# Example:
#   scripts/testing/run_and_log.sh A1 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 examples/magnetostatics/01_straight_wire.py'"

DRY_RUN="${FEM_SOLVER_DRY_RUN:-0}"

if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=1
  shift
fi

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 [--dry-run] <chunk_id> <command...>"
  echo "       FEM_SOLVER_DRY_RUN=1 $0 <chunk_id> <command...>"
  exit 2
fi

CHUNK_ID="$1"
shift
CMD="$*"

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
LOG_DIR="$ROOT_DIR/docs/testing/logs"
INDEX_FILE="$ROOT_DIR/docs/testing/test-results.md"
mkdir -p "$LOG_DIR"

# Ensure docker compose can always find project config for this repo layout.
DEFAULT_COMPOSE_FILE="$ROOT_DIR/docker/docker-compose.yml"
if [[ -z "${COMPOSE_FILE:-}" && -f "$DEFAULT_COMPOSE_FILE" ]]; then
  export COMPOSE_FILE="$DEFAULT_COMPOSE_FILE"
fi

ensure_index_file() {
  if [[ ! -f "$INDEX_FILE" ]]; then
    cat > "$INDEX_FILE" <<'EOF'
# FEM-EM Manual Test Results

This file is append-only and records test runs executed by the human operator.
Each run has a full log in `docs/testing/logs/`.

| UTC Timestamp | Chunk | Commit | Elapsed (s) | Env | Exit | Log |
|---|---|---|---:|---|---:|---|
EOF
    return
  fi

  python3 - "$INDEX_FILE" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
text = path.read_text()
legacy_header = "| UTC Timestamp | Chunk | Exit | Log |"
legacy_sep = "|---|---|---:|---|"
new_header = "| UTC Timestamp | Chunk | Commit | Elapsed (s) | Env | Exit | Log |"
new_sep = "|---|---|---|---:|---|---:|---|"

if new_header in text:
    sys.exit(0)

if legacy_header not in text:
    sys.exit(0)

lines = text.splitlines()
out = []
i = 0
while i < len(lines):
    line = lines[i]
    if line == legacy_header:
        out.append(new_header)
        if i + 1 < len(lines) and lines[i + 1] == legacy_sep:
            out.append(new_sep)
            i += 2
        else:
            out.append(new_sep)
            i += 1
        while i < len(lines):
            row = lines[i]
            if not row.startswith("|"):
                out.append(row)
                i += 1
                continue
            parts = [p.strip() for p in row.strip().split("|")[1:-1]]
            if len(parts) == 4:
                utc_ts, chunk, exit_code, log_ref = parts
                out.append(f"| {utc_ts} | {chunk} | `unknown` | unknown | `legacy` | {exit_code} | {log_ref} |")
            else:
                out.append(row)
            i += 1
        break
    out.append(line)
    i += 1

path.write_text("\n".join(out) + "\n")
PY
}

TS="$(date -u +%Y%m%dT%H%M%SZ)"
SAFE_CHUNK="${CHUNK_ID//[^a-zA-Z0-9._-]/_}"
LOG_FILE="$LOG_DIR/${TS}_${SAFE_CHUNK}.log"

GIT_COMMIT_FULL="$(git -C "$ROOT_DIR" rev-parse HEAD 2>/dev/null || echo unknown)"
DOCKER_COMPOSE_VERSION="$(docker compose version 2>/dev/null | head -n 1 || echo unavailable)"
OS_META="$(uname -srm 2>/dev/null || echo unknown)"
PYTHON_META="$(python3 --version 2>/dev/null || echo unavailable)"
HOST_META="$(hostname 2>/dev/null || echo unknown)"
COMPOSE_META="${COMPOSE_FILE:-<unset>}"
ENV_META="host=${HOST_META};os=${OS_META};python=${PYTHON_META};compose=$(basename "$COMPOSE_META")"
ENV_META="${ENV_META//|//}"

# Dry-run mode: just print what would be executed, don't run anything
if [[ "$DRY_RUN" == "1" ]]; then
  echo "[DRY RUN] Would execute test for chunk: $CHUNK_ID"
  echo "[DRY RUN] Command: $CMD"
  echo "[DRY RUN] Commit: $GIT_COMMIT_FULL"
  echo "[DRY RUN] Env: $ENV_META"
  echo "[DRY RUN] Log would be written to: $LOG_FILE"
  echo ""
  echo "To actually run this test, execute:"
  echo "  $0 $CHUNK_ID \"$CMD\""
  exit 0
fi

START_EPOCH="$(date -u +%s)"

{
  echo "# Test Run"
  echo "- Chunk: $CHUNK_ID"
  echo "- Commit: $GIT_COMMIT_FULL"
  echo "- Time (UTC): $(date -u '+%Y-%m-%d %H:%M:%S')"
  echo "- Host: $HOST_META"
  echo "- OS: $OS_META"
  echo "- Python: $PYTHON_META"
  echo "- Docker Compose: $DOCKER_COMPOSE_VERSION"
  echo "- Root Dir: $ROOT_DIR"
  echo "- Working Dir: $(pwd)"
  echo "- COMPOSE_FILE: $COMPOSE_META"
  echo "- Command: $CMD"
  echo ""
  echo "## Preflight"
  echo '\$ pwd'
  pwd
  echo ""
  echo '\$ ls -la docker'
  ls -la "$ROOT_DIR/docker" 2>&1 || true
  echo ""
  echo '\$ docker compose version'
  docker compose version 2>&1 || true
  echo ""
  echo '\$ docker compose config --services'
  docker compose config --services 2>&1 || true
  echo ""
  echo "## Output"
} > "$LOG_FILE"

set +e
(
  cd "$ROOT_DIR"
  set -x
  bash -lc "$CMD"
) >> "$LOG_FILE" 2>&1
STATUS=$?
set -e

END_EPOCH="$(date -u +%s)"
ELAPSED_SECONDS=$((END_EPOCH - START_EPOCH))

{
  echo ""
  echo "## Exit"
  echo "- Status: $STATUS"
  echo "- Elapsed (s): $ELAPSED_SECONDS"
} >> "$LOG_FILE"

ensure_index_file
printf '| %s | %s | `%s` | %s | `%s` | %s | `%s` |\n' \
  "$(date -u '+%Y-%m-%d %H:%M:%S')" \
  "$CHUNK_ID" \
  "$GIT_COMMIT_FULL" \
  "$ELAPSED_SECONDS" \
  "$ENV_META" \
  "$STATUS" \
  "$(basename "$LOG_FILE")" >> "$INDEX_FILE"

echo "Log written: $LOG_FILE"
echo "Index updated: $INDEX_FILE"
exit "$STATUS"
