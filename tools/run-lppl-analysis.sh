#!/bin/bash
# Run full LPPL analysis on all Olympus scan logs, update summary.
# Usage: bash tools/run-lppl-analysis.sh
set -e

echo "Syncing lppl-v3 to Olympus..."
scp tools/scan-lppl-v3.py olympus:/tmp/scan-lppl-v3.py

echo "Running LPPL v3 on all logs..."
ssh olympus "python3 /tmp/scan-lppl-v3.py --all-logs /tmp --json /tmp/lppl_v3_latest.json"

echo "Fetching results..."
scp olympus:/tmp/lppl_v3_latest.json /tmp/lppl_v3_latest.json

echo "Running convergence analysis..."
python3 tools/scan-lppl-convergence.py /tmp/lppl_v3_latest.json

echo "Done."
