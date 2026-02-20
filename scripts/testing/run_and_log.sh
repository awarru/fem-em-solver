#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   scripts/testing/run_and_log.sh <chunk_id> <command...>
# Example:
#   scripts/testing/run_and_log.sh A1 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 examples/magnetostatics/01_straight_wire.py'"

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <chunk_id> <command...>"
  exit 2
fi

CHUNK_ID="$1"
shift
CMD="$*"

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
LOG_DIR="$ROOT_DIR/docs/testing/logs"
INDEX_FILE="$ROOT_DIR/docs/testing/test-results.md"
mkdir -p "$LOG_DIR"

TS="$(date -u +%Y%m%dT%H%M%SZ)"
SAFE_CHUNK="${CHUNK_ID//[^a-zA-Z0-9._-]/_}"
LOG_FILE="$LOG_DIR/${TS}_${SAFE_CHUNK}.log"

{
  echo "# Test Run"
  echo "- Chunk: $CHUNK_ID"
  echo "- Time (UTC): $(date -u '+%Y-%m-%d %H:%M:%S')"
  echo "- Host: $(hostname)"
  echo "- Command: $CMD"
  echo ""
  echo "## Output"
} > "$LOG_FILE"

set +e
bash -lc "$CMD" >> "$LOG_FILE" 2>&1
STATUS=$?
set -e

{
  echo ""
  echo "## Exit"
  echo "- Status: $STATUS"
} >> "$LOG_FILE"

if [[ ! -f "$INDEX_FILE" ]]; then
  cat > "$INDEX_FILE" <<'EOF'
# FEM-EM Manual Test Results

This file is append-only and records test runs executed by the human operator.
Each run has a full log in `docs/testing/logs/`.

| UTC Timestamp | Chunk | Exit | Log |
|---|---|---:|---|
EOF
fi

printf '| %s | %s | %s | `%s` |\n' "$(date -u '+%Y-%m-%d %H:%M:%S')" "$CHUNK_ID" "$STATUS" "$(basename "$LOG_FILE")" >> "$INDEX_FILE"

echo "Log written: $LOG_FILE"
echo "Index updated: $INDEX_FILE"
exit "$STATUS"
