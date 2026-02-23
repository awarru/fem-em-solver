#!/usr/bin/env bash
set -euo pipefail

# Run one or more example scripts inside docker compose without entering container.
#
# Usage:
#   ./scripts/run_examples.sh --list
#   ./scripts/run_examples.sh -e 1
#   ./scripts/run_examples.sh -e 1,3 -n 4
#   ./scripts/run_examples.sh -e all -n 2
#
# Options:
#   -e, --example   Example selection by number (e.g. 1 or 01), CSV (1,3), or 'all'
#   -n, --nproc     MPI process count (default: 2)
#   --list          List available examples and exit

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
EXAMPLES_DIR="$ROOT_DIR/examples/magnetostatics"
DEFAULT_COMPOSE_FILE="$ROOT_DIR/docker/docker-compose.yml"

if [[ -z "${COMPOSE_FILE:-}" && -f "$DEFAULT_COMPOSE_FILE" ]]; then
  export COMPOSE_FILE="$DEFAULT_COMPOSE_FILE"
fi

NPROC=2
EXAMPLE_SPEC=""
MODE="run"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -e|--example)
      EXAMPLE_SPEC="${2:-}"
      [[ -n "$EXAMPLE_SPEC" ]] || { echo "--example requires a value" >&2; exit 2; }
      shift 2
      ;;
    -n|--nproc)
      NPROC="${2:-}"
      [[ "$NPROC" =~ ^[1-9][0-9]*$ ]] || { echo "--nproc must be a positive integer" >&2; exit 2; }
      shift 2
      ;;
    --list)
      MODE="list"
      shift
      ;;
    -h|--help)
      sed -n '1,35p' "$0"
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

mapfile -t AVAILABLE < <(find "$EXAMPLES_DIR" -maxdepth 1 -type f -name '*.py' | sort)

if [[ ${#AVAILABLE[@]} -eq 0 ]]; then
  echo "No examples found in $EXAMPLES_DIR" >&2
  exit 1
fi

if [[ "$MODE" == "list" ]]; then
  echo "Available examples:"
  for p in "${AVAILABLE[@]}"; do
    b="$(basename "$p")"
    num="${b%%_*}"
    # normalize display: 01 -> 1
    num_display="$((10#$num))"
    echo "  $num_display -> ${p#$ROOT_DIR/}"
  done
  exit 0
fi

if [[ -z "$EXAMPLE_SPEC" ]]; then
  echo "Missing --example selection. Use --list to see available options." >&2
  exit 2
fi

SELECTED=()
if [[ "$EXAMPLE_SPEC" == "all" ]]; then
  SELECTED=("${AVAILABLE[@]}")
else
  IFS=',' read -r -a TOKENS <<< "$EXAMPLE_SPEC"
  declare -A SEEN=()

  for t in "${TOKENS[@]}"; do
    token="$(echo "$t" | xargs)"
    [[ "$token" =~ ^[0-9]+$ ]] || { echo "Invalid example token: '$token'" >&2; exit 2; }

    num_padded=$(printf "%02d" "$((10#$token))")
    matches=("$EXAMPLES_DIR/${num_padded}"_*.py)

    if [[ ! -e "${matches[0]}" ]]; then
      echo "No example found for number: $token" >&2
      exit 2
    fi

    # if multiple match weirdly, run all but de-dup
    for m in "${matches[@]}"; do
      [[ -f "$m" ]] || continue
      if [[ -z "${SEEN[$m]:-}" ]]; then
        SELECTED+=("$m")
        SEEN[$m]=1
      fi
    done
  done
fi

cd "$ROOT_DIR"

echo "Running ${#SELECTED[@]} example(s) with mpiexec -n $NPROC"
echo "COMPOSE_FILE=${COMPOSE_FILE:-<unset>}"

for ex in "${SELECTED[@]}"; do
  rel="${ex#$ROOT_DIR/}"
  echo
  echo "==> $rel"
  docker compose exec fem-em-solver bash -lc "cd /workspace && PYTHONPATH=/workspace/src mpiexec -n $NPROC python3 $rel"
done

echo
echo "Done."
