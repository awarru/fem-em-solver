#!/usr/bin/env bash
set -euo pipefail

# Project-root test entrypoint.
# Mirrors scripts/run_tests.sh so you can run from repo root with a single command.

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
exec "$ROOT_DIR/scripts/testing/run_pending_tests.sh" "$@"
