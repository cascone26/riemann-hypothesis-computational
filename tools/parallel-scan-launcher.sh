#!/bin/bash
# parallel-scan-launcher.sh
# Launch multiple dense_sys scans in parallel on Olympus,
# each covering a different lsq range.
# Usage: ./parallel-scan-launcher.sh <base_lsq> <n_parallel>
# e.g.: ./parallel-scan-launcher.sh 28000 4
# Launches ranges: [28001-28500], [30001-30500], [32001-32500], [34001-34500]

OLYMPUS="olympus"
BASE_LSQ=${1:-28000}
N_PARALLEL=${2:-4}
STEP=500
GAP=2000  # lsq gap between parallel ranges

echo "Launching $N_PARALLEL parallel scans from lsq $BASE_LSQ"
echo "Ranges:"

for i in $(seq 0 $((N_PARALLEL-1))); do
    START=$((BASE_LSQ + 1 + i*GAP))
    END=$((START + STEP - 1))
    LABEL="PAR${i}${BASE_LSQ}"
    echo "  Range $i: lsq $START - $END (label $LABEL)"
done

echo ""
echo "Creating scripts on Olympus..."

for i in $(seq 0 $((N_PARALLEL-1))); do
    START=$((BASE_LSQ + 1 + i*GAP))
    END=$((START + STEP))
    LABEL="P${i}B"
    SCRIPT="/tmp/dense_sys_${LABEL}.py"
    LOGFILE="/tmp/dense_sys_${LABEL}.log"

    ssh $OLYMPUS "python3 - <<'PYEOF'
import re
content = open('/tmp/dense_sys_BC.py').read()
# Replace range
m = re.search(r'range\((\d+),\s*(\d+)\)', content)
if m:
    content = content.replace(m.group(0), f'range($START, $END)')
# Replace log label in output (replace DONE BC with DONE $LABEL)
content = content.replace('DONE BC', 'DONE $LABEL')
open('$SCRIPT', 'w').write(content)
print(f'Created $SCRIPT: range($START, $END)')
PYEOF"

    ssh $OLYMPUS "nohup python3 $SCRIPT > $LOGFILE 2>&1 &"
    echo "  Launched $LABEL (range $START-$END)"
done

echo ""
echo "All launched. Monitor with:"
echo "  ssh $OLYMPUS 'tail -1 /tmp/dense_sys_P*B.log'"
