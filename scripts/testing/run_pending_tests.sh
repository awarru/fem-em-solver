#!/usr/bin/env bash
set -euo pipefail

# Run manual ROADMAP pending tests by executing commands listed in
# docs/testing/pending-tests.md (each wrapped with scripts/testing/run_and_log.sh).
#
# Usage:
#   ./scripts/testing/run_pending_tests.sh                 # run all pending manual tests
#   ./scripts/testing/run_pending_tests.sh --list          # list discovered chunk commands + recommendation order
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
RESULTS_FILE="$ROOT_DIR/docs/testing/test-results.md"

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

mapfile -t RAW_COMMAND_LINES < <(
  grep -E 'scripts/testing/run_and_log\.sh[[:space:]]+[A-Za-z0-9._-]+[[:space:]]+".*"' "$PENDING_FILE" \
    | sed -E 's/^.*(scripts\/testing\/run_and_log\.sh[[:space:]]+[A-Za-z0-9._-]+[[:space:]]+".*")$/\1/'
)

if [[ ${#RAW_COMMAND_LINES[@]} -eq 0 ]]; then
  echo "No pending manual test commands found in $PENDING_FILE"
  exit 0
fi

declare -A COMMAND_BY_CHUNK=()
for cmd in "${RAW_COMMAND_LINES[@]}"; do
  chunk="$(awk '{print $2}' <<<"$cmd")"
  if [[ -n "$chunk" ]]; then
    # Keep the last command for a chunk (newest entry wins).
    COMMAND_BY_CHUNK["$chunk"]="$cmd"
  fi
done

declare -A STATUS_BY_CHUNK=()
while IFS='|' read -r chunk status; do
  [[ -z "$chunk" ]] && continue
  STATUS_BY_CHUNK["$chunk"]="$status"
done < <(
  awk -F'|' '
    /^\|/ {
      chunk=$2; status=$3
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", chunk)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", status)
      if (chunk ~ /^[A-Za-z][0-9]+$/) {
        print chunk "|" status
      }
    }
  ' "$PENDING_FILE"
)

declare -A EXIT_BY_CHUNK=()
if [[ -f "$RESULTS_FILE" ]]; then
  while IFS='|' read -r chunk exit_code; do
    [[ -z "$chunk" ]] && continue
    EXIT_BY_CHUNK["$chunk"]="$exit_code"
  done < <(
    awk -F'|' '
      /^\|[[:space:]]*[0-9]{4}-[0-9]{2}-[0-9]{2}/ {
        chunk=$3; exit_code=$7
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", chunk)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", exit_code)
        if (chunk ~ /^[A-Za-z][0-9]+$/) {
          latest[chunk]=exit_code
        }
      }
      END {
        for (k in latest) {
          print k "|" latest[k]
        }
      }
    ' "$RESULTS_FILE"
  )
fi

mapfile -t ALL_CHUNKS < <(printf '%s\n' "${!COMMAND_BY_CHUNK[@]}" | sort -V)

if [[ -n "$CHUNK_FILTER" ]]; then
  if [[ -z "${COMMAND_BY_CHUNK[$CHUNK_FILTER]+x}" ]]; then
    echo "No pending manual test command found for chunk: $CHUNK_FILTER" >&2
    exit 2
  fi
  ALL_CHUNKS=("$CHUNK_FILTER")
fi

mapfile -t COMMAND_LINES < <(
  for chunk in "${ALL_CHUNKS[@]}"; do
    echo "${COMMAND_BY_CHUNK[$chunk]}"
  done
)

build_recommendation_table() {
  local chunk status last_exit score reason

  for chunk in "${ALL_CHUNKS[@]}"; do
    status="${STATUS_BY_CHUNK[$chunk]:-unknown}"
    last_exit="${EXIT_BY_CHUNK[$chunk]:-unknown}"

    score=300
    reason="queued"

    if [[ "$status" == *"BLOCKED"* ]]; then
      score=100
      reason="unblocker"
    elif [[ "$status" == *"AWAITING-HUMAN-TEST"* ]]; then
      score=200
      reason="awaiting human verification"
    elif [[ "$status" == *"COMPLETE"* ]]; then
      score=500
      reason="already complete"
    fi

    if [[ "$last_exit" =~ ^[0-9]+$ ]] && [[ "$last_exit" -ne 0 ]]; then
      score=$((score - 20))
      reason="$reason; last run exited $last_exit"
    fi

    # Stable tie-break by chunk id.
    printf '%04d|%s|%s|%s|%s\n' "$score" "$chunk" "$status" "$last_exit" "$reason"
  done | sort -t'|' -k1,1n -k2,2V
}

if [[ "$MODE" == "list" ]]; then
  echo "Discovered manual test commands (latest per chunk):"
  for cmd in "${COMMAND_LINES[@]}"; do
    echo "  $cmd"
  done

  echo ""
  echo "Recommended next test order:"
  rank=1
  while IFS='|' read -r _score chunk status last_exit reason; do
    cmd="${COMMAND_BY_CHUNK[$chunk]}"
    echo "  $rank) $chunk [$status]"
    echo "     reason: $reason"
    echo "     last exit: $last_exit"
    echo "     command: $cmd"
    rank=$((rank + 1))
  done < <(build_recommendation_table)

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
  echo "To see the list of pending commands and recommended order:"
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
