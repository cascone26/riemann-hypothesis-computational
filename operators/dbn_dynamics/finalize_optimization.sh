#!/bin/bash
# Finalize optimization: validate results and deploy

cd /Users/scones/projects/riemann/operators/dbn_dynamics

PYTHON=/Users/scones/projects/riemann/flint_venv/bin/python3
FINAL="arb_certify2_final.py"
MAIN="arb_certify2.py"

echo "=== OPTIMIZATION FINALIZATION ==="
echo

# Wait for the final run to complete (called externally by monitor)
# This script is run AFTER the monitor reports completion

if [ ! -f "${FINAL}" ]; then
    echo "ERROR: ${FINAL} not found"
    exit 1
fi

echo "Step 1: Validating optimized version..."
if grep -q "*** CERTIFIED PD" /private/tmp/claude-501/-Users-scones/76e94d23-a249-471d-9af6-141ebd7506c0/tasks/b0wfj93fj.output 2>/dev/null; then
    echo "✓ CERTIFICATION PASSED"
else
    echo "✗ Certification check failed - checking file for errors"
    cat /private/tmp/claude-501/-Users-scones/76e94d23-a249-471d-9af6-141ebd7506c0/tasks/b0wfj93fj.output 2>/dev/null | tail -30
    exit 1
fi

echo
echo "Step 2: Extracting results..."
OUTPUT=$(/private/tmp/claude-501/-Users-scones/76e94d23-a249-471d-9af6-141ebd7506c0/tasks/b0wfj93fj.output 2>/dev/null)
RUNTIME=$(echo "$OUTPUT" | grep "runtime=" | sed 's/.*runtime=\([0-9.]*\)s.*/\1/')
echo "Runtime: ${RUNTIME}s"

if [ -z "$RUNTIME" ]; then
    echo "Could not extract runtime, checking output:"
    tail -15 /private/tmp/claude-501/-Users-scones/76e94d23-a249-471d-9af6-141ebd7506c0/tasks/b0wfj93fj.output 2>/dev/null
else
    echo "✓ Runtime extracted: ${RUNTIME}s"
fi

echo
echo "Step 3: Deploying optimized version..."
cp "${FINAL}" "${MAIN}"
echo "✓ Copied ${FINAL} -> ${MAIN}"

echo
echo "=== FINAL OPTIMIZED RESULTS ==="
echo
$PYTHON "${MAIN}"
