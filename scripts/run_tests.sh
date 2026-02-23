#!/usr/bin/env bash
set -euo pipefail

# Convenience wrapper: run pending human-manual tests from docs/testing/pending-tests.md
# with automatic logging to docs/testing/logs and docs/testing/test-results.md.
#
# Usage:
#   ./scripts/run_tests.sh
#   ./scripts/run_tests.sh --list
#   ./scripts/run_tests.sh --chunk E3

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
exec "$ROOT_DIR/scripts/testing/run_pending_tests.sh" "$@"
