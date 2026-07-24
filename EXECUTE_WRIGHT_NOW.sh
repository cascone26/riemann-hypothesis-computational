#!/bin/bash
# Execute full Wright cycle on all 4 sorries in RobinInequality.lean
# Call this the MOMENT the lake build completes

set -e
cd /Users/scones/projects/riemann

echo ""
echo "================================================================================"
echo " WRIGHT HARNESS - FULL EXECUTION CYCLE"
echo " Time: $(date)"
echo "================================================================================"
echo ""

# Quick verification that build succeeded
echo "Verifying lake build status..."
cd RiemannLean
if lake build 2>&1 | tail -5 | grep -q "error\|Error"; then
  echo "ERROR: Build has failures"
  exit 1
fi
echo "✓ Build successful"
cd ..

# Run the Wright test
echo ""
echo "Executing Wright cycle on RobinInequality.lean..."
echo "--------"
python3 bot/wright_lean_harness.py 2>&1 | tee wright_execution.log

echo ""
echo "================================================================================"
echo " Wright cycle execution complete"
echo " Output saved to: wright_execution.log"
echo "================================================================================"
