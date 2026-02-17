#!/usr/bin/env bash
# Build documentation

cd "$(dirname "$0")/.."

echo "Building documentation..."
mkdocs build

echo ""
echo "Documentation built in site/"
echo "Run 'mkdocs serve' to preview locally"
