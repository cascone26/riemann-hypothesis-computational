#!/bin/bash
# jump-scan.sh — launch scans at a fresh lsq range to escape blob deserts.
# Usage: ./jump-scan.sh <start_lsq> [n_ranges=4]
# e.g.: ./jump-scan.sh 35000 4
# Creates and launches 4 parallel 500-step scans at well-separated starting points.

OLYMPUS="olympus"
START=${1:-35000}
N=${2:-4}
STEP=500
GAP=5000  # Spread out ranges to sample different lsq territory

echo "Launching $N jump scans starting at lsq $START"
echo "Spacing: $GAP lsq between each"
echo ""

# Need a template script — use BC as template (last known working script)
TEMPLATE_SCRIPT="/tmp/dense_sys_BC.py"

for i in $(seq 0 $((N-1))); do
    RANGE_START=$((START + i*GAP + 1))
    RANGE_END=$((RANGE_START + STEP))
    LABEL="J${i}C"
    SCRIPT="/tmp/dense_sys_${LABEL}.py"
    LOG="/tmp/dense_sys_${LABEL}.log"

    echo "  J${i}: lsq $RANGE_START-$((RANGE_END-1)) -> $SCRIPT"

    ssh $OLYMPUS "python3 -c \"
import re
content = open('$TEMPLATE_SCRIPT').read()
m = re.search(r'range\((\d+),\s*(\d+)\)', content)
if m:
    content = content.replace(m.group(0), 'range($RANGE_START, $RANGE_END)')
# Fix DONE label
content = content.replace('DONE BC', 'DONE $LABEL')
open('$SCRIPT', 'w').write(content)
print('Created: range($RANGE_START, $RANGE_END)')
\"" && ssh $OLYMPUS "nohup python3 $SCRIPT > $LOG 2>&1 &" && echo "    Launched." || echo "    FAILED."
done

echo ""
echo "Monitor: ssh $OLYMPUS \"for f in /tmp/dense_sys_J*.log; do echo \$f: \$(tail -1 \$f 2>/dev/null); done\""
