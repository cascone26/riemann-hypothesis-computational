#!/bin/bash
# Portable never-stop wrapper. Halts ONLY on PROOF_FOUND.marker. Heartbeat for liveness (no phone spam).
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$DIR")"
cd "$DIR"
while [ ! -f "$ROOT/results/PROOF_FOUND.marker" ]; do
  date +%s > "$DIR/heartbeat"                 # liveness (checkable, silent)
  python3 rh_engine.py 1 >> engine_run.log 2>&1
  sleep 60
done
echo "PROOF_FOUND — machine halted $(date)" >> engine_run.log
