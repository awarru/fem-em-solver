#!/usr/bin/env bash
set -euo pipefail

# Run manual ROADMAP pending tests by executing commands listed in
# docs/testing/pending-tests.md (each wrapped with scripts/testing/run_and_log.sh).
#
# Usage:
#   ./scripts/testing/run_pending_tests.sh            # run all pending manual tests
#   ./scripts/testing/run_pending_tests.sh --list     # list discovered chunk commands
#   ./scripts/testing/run_pending_tests.sh --chunk E3 # run only one chunk

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
PENDING_FILE="$ROOT_DIR/docs/testing/pending-tests.md"

if [[ ! -f "$PENDING_FILE" ]]; then
  echo "Missing pending test file: $PENDING_FILE" >&2
  exit 2
fi

MODE="run"
CHUNK_FILTER=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --list)
      MODE="list"
      shift
      ;;
    --chunk)
      CHUNK_FILTER="${2:-}"
      if [[ -z "$CHUNK_FILTER" ]]; then
        echo "--chunk requires a chunk id (e.g. E3)" >&2
        exit 2
      fi
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: $0 [--list] [--chunk <ID>]" >&2
      exit 2
      ;;
  esac
done

mapfile -t COMMAND_LINES < <(
  grep -E '^[[:space:]]*scripts/testing/run_and_log\.sh[[:space:]]+[A-Za-z0-9._-]+[[:space:]]+".*"[[:space:]]*$' "$PENDING_FILE" \
    | sed -E 's/^[[:space:]]*//'
)

if [[ ${#COMMAND_LINES[@]} -eq 0 ]]; then
  echo "No pending manual test commands found in $PENDING_FILE"
  exit 0
fi

if [[ -n "$CHUNK_FILTER" ]]; then
  mapfile -t COMMAND_LINES < <(
    printf '%s\n' "${COMMAND_LINES[@]}" | awk -v c="$CHUNK_FILTER" '$2 == c { print }'
  )

  if [[ ${#COMMAND_LINES[@]} -eq 0 ]]; then
    echo "No pending manual test command found for chunk: $CHUNK_FILTER" >&2
    exit 2
  fi
fi

if [[ "$MODE" == "list" ]]; then
  printf '%s\n' "${COMMAND_LINES[@]}"
  exit 0
fi

cd "$ROOT_DIR"

for cmd in "${COMMAND_LINES[@]}"; do
  echo
  echo "==> Running: $cmd"
  bash -lc "$cmd"
done

echo
echo "Done. Logs: $ROOT_DIR/docs/testing/logs"
echo "Index: $ROOT_DIR/docs/testing/test-results.md"
