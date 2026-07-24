#!/bin/bash
# Complete Wright harness test - runs when lake build is ready

set -e

cd /Users/scones/projects/riemann

echo ""
echo "================================================================================"
echo "WRIGHT HARNESS TEST — REAL BUILD CYCLE"
echo "================================================================================"
echo "Test execution time: $(date)"
echo ""

# Verify lake build succeeded
if ! ~/.elan/bin/lake build 2>&1 | tail -5; then
  echo "ERROR: lake build failed!"
  exit 1
fi

echo ""
echo "✓ Lake build successful"
echo ""
echo "Running Wright harness on RobinInequality.lean sorries..."
echo ""

# Run the harness
if python3 bot/wright_lean_harness.py; then
  echo ""
  echo "✓ Wright harness test completed successfully"
else
  echo ""
  echo "✗ Wright harness encountered errors"
  exit 1
fi

echo ""
echo "================================================================================"
echo "Test execution complete at: $(date)"
echo "================================================================================"
