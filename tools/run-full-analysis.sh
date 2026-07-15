#!/bin/bash
# Full analysis pipeline — run on Olympus after each scan batch completes.
# Usage: ./run-full-analysis.sh
# 1. Regenerates lppl_results_all.json from all dense_sys_*.log files
# 2. Runs omega histogram + clustering test
# 3. Runs phase coherence test
# 4. Runs m-exponent test
# 5. Shows convergence summary

set -e
OLYMPUS="olympus"

echo "=== Step 1: Regenerate all-logs JSON ==="
ssh $OLYMPUS "python3 /tmp/lppl_v3.py --all-logs /tmp --json /tmp/lppl_results_all.json"

echo ""
echo "=== Step 2: Pull results to Mac ==="
scp $OLYMPUS:/tmp/lppl_results_all.json /tmp/lppl_results_all_local.json

echo ""
echo "=== Step 3: Omega histogram + clustering ==="
ssh $OLYMPUS "python3 /tmp/omega_histogram.py"

echo ""
echo "=== Step 4: Phase coherence ==="
ssh $OLYMPUS "python3 /tmp/phase_test.py"

echo ""
echo "=== Step 5: m-exponent test ==="
ssh $OLYMPUS "python3 /tmp/m_test.py"

echo ""
echo "=== Step 6: Summary stats ==="
ssh $OLYMPUS "python3 /tmp/lppl_stat2.py"

echo ""
echo "=== Scan status ==="
ssh $OLYMPUS "ls -la /tmp/dense_sys_B*.log | tail -5"
ssh $OLYMPUS "tail -1 /tmp/dense_sys_BC.log 2>/dev/null || echo 'BC not running'"
