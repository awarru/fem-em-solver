#!/usr/bin/env bash
# Run all tests

cd "$(dirname "$0")/.."

echo "Running unit tests..."
pytest tests/unit -v -x

echo ""
echo "Running integration tests..."
pytest tests/integration -v -x

echo ""
echo "Running validation tests..."
pytest tests/validation -v -x

echo ""
echo "All tests completed!"
