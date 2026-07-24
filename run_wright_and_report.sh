#!/bin/bash
# Run Wright harness and capture real build output

set -e

cd /Users/scones/projects/riemann

echo "================================================================================"
echo "WRIGHT HARNESS EXECUTION — REAL BUILD CYCLE"
echo "================================================================================"
echo "Time: $(date)"
echo ""

# Run the test harness
echo "Running: python3 test_wright_cycle.py"
echo "--------"
python3 test_wright_cycle.py 2>&1 | tee wright_cycle_output.log

echo ""
echo "================================================================================"
echo "Cycle complete. Output saved to: wright_cycle_output.log"
echo "================================================================================"
