#!/usr/bin/env bash
set -euo pipefail

# Run manual ROADMAP pending tests by executing commands listed in
# docs/testing/pending-tests.md (each wrapped with scripts/testing/run_and_log.sh).
#
# Usage:
#   ./scripts/testing/run_pending_tests.sh                 # run all pending manual tests
#   ./scripts/testing/run_pending_tests.sh --list          # list discovered chunk commands
#   ./scripts/testing/run_pending_tests.sh --chunk E3      # run only one chunk
#   ./scripts/testing/run_pending_tests.sh --smoke         # run lightweight smoke matrix (cron-safe)
#   ./scripts/testing/run_pending_tests.sh --dry-run       # show what would be run (cron-safe)
#   FEM_SOLVER_DRY_RUN=1 ./scripts/testing/run_pending_tests.sh
#
# Cron/Automation Note:
#   This script defaults to EXECUTING tests. For automated environments (like cron),
#   use --dry-run or set FEM_SOLVER_DRY_RUN=1 to prevent accidental CPU-heavy execution.

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
PENDING_FILE="$ROOT_DIR/docs/testing/pending-tests.md"

SMOKE_TEST_TARGETS=(
  "tests/unit/test_analytical_lightweight.py"
  "tests/solver/test_tolerance_policy.py"
  "tests/validation/test_tolerance_policy.py"
  "tests/ports/test_port_definition.py"
)

if [[ ! -f "$PENDING_FILE" ]]; then
  echo "Missing pending test file: $PENDING_FILE" >&2
  exit 2
fi

MODE="run"
CHUNK_FILTER=""
DRY_RUN="${FEM_SOLVER_DRY_RUN:-0}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --list)
      MODE="list"
      shift
      ;;
    --smoke)
      MODE="smoke"
      shift
      ;;
    --dry-run)
      DRY_RUN=1
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
      echo "Usage: $0 [--list] [--smoke] [--dry-run] [--chunk <ID>]" >&2
      exit 2
      ;;
  esac
done

if [[ "$MODE" == "smoke" && -n "$CHUNK_FILTER" ]]; then
  echo "--smoke cannot be combined with --chunk" >&2
  exit 2
fi

if [[ "$MODE" == "smoke" ]]; then
  if [[ "$DRY_RUN" == "1" ]]; then
    echo "=== DRY RUN MODE (no tests will be executed) ==="
    echo ""
    echo "Smoke matrix command:"
    echo "  python3 -m pytest ${SMOKE_TEST_TARGETS[*]} -v"
    exit 0
  fi

  cd "$ROOT_DIR"
  echo "Running lightweight smoke matrix (no heavy FEM/mesh/solver cases):"
  for target in "${SMOKE_TEST_TARGETS[@]}"; do
    echo "  - $target"
  done
  python3 -m pytest "${SMOKE_TEST_TARGETS[@]}" -v
  exit 0
fi

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

# Dry-run mode: show what would be executed without running
if [[ "$DRY_RUN" == "1" ]]; then
  echo "=== DRY RUN MODE (no tests will be executed) ==="
  echo ""
  echo "The following tests are pending and need to be run manually:"
  echo ""
  for cmd in "${COMMAND_LINES[@]}"; do
    echo "  $cmd"
  done
  echo ""
  echo "To run these tests manually:"
  echo "  $0"
  echo ""
  echo "To run a specific chunk:"
  echo "  $0 --chunk <ID>"
  echo ""
  echo "To run the lightweight smoke matrix:"
  echo "  $0 --smoke"
  echo ""
  echo "To see the list of pending commands:"
  echo "  $0 --list"
  exit 0
fi

"$ROOT_DIR/scripts/testing/preflight.sh"

cd "$ROOT_DIR"

for cmd in "${COMMAND_LINES[@]}"; do
  echo
  echo "==> Running: $cmd"
  bash -lc "$cmd"
done

echo
echo "Done. Logs: $ROOT_DIR/docs/testing/logs"
echo "Index: $ROOT_DIR/docs/testing/test-results.md"
