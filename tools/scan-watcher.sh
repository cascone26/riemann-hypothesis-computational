#!/usr/bin/env bash
# scan-watcher.sh — Monitor a dense_sys scan log and send ntfy on key events
# Usage: ./scan-watcher.sh /tmp/dense_sys_AS.log [ntfy-topic]
# Runs on Olympus. Zero Claude tokens.

LOG="${1:-/tmp/dense_sys_AS.log}"
TOPIC="${2:-jacob-cobo-80c6d2b9e9c4}"
NTFY="https://ntfy.sh/$TOPIC"

BC_COUNT=0
ZONE_STATE="unknown"  # "blob" or "zone"
PREV_DEEP=""

notify() {
    curl -s -d "$1" "$NTFY" -H "Title: Riemann Scan" -H "Priority: high" > /dev/null 2>&1
}

echo "Watching $LOG for events → ntfy $TOPIC"
notify "Scan watcher started: $LOG"

tail -n +1 -F "$LOG" | while IFS= read -r line; do
    # Scan complete
    if echo "$line" | grep -qE "^DONE "; then
        SCAN_NAME=$(echo "$line" | awk '{print $2}')
        notify "DONE: $SCAN_NAME scan complete"
        echo "DONE: $SCAN_NAME"
        break
    fi

    # BC event (RECORD-CLASS)
    if echo "$line" | grep -q "RECORD-CLASS"; then
        LSQ=$(echo "$line" | awk '{print $1}')
        L1=$(echo "$line" | awk '{print $2}')
        L2=$(echo "$line" | awk '{print $3}')
        BC_COUNT=$((BC_COUNT+1))
        notify "BC event lsq $LSQ: λ₁=$L1 λ₂=$L2 (total BCs this watch: $BC_COUNT)"
        echo "BC at lsq $LSQ"
        continue
    fi

    # Zone entry (DEEP line after non-DEEP)
    if echo "$line" | grep -q "<--DEEP"; then
        if [ "$ZONE_STATE" = "blob" ]; then
            LSQ=$(echo "$line" | awk '{print $1}')
            L1=$(echo "$line" | awk '{print $2}')
            notify "ZONE ENTRY lsq $LSQ: λ₁=$L1"
            echo "Zone entry at lsq $LSQ"
        fi
        ZONE_STATE="zone"
        continue
    fi

    # Blob entry (non-DEEP line after DEEP)
    if [ "$ZONE_STATE" = "zone" ] && ! echo "$line" | grep -q "<--"; then
        LSQ=$(echo "$line" | awk '{print $1}')
        L1=$(echo "$line" | awk '{print $2}')
        if [ -n "$LSQ" ] && echo "$L1" | grep -qE '^-[0-9]'; then
            notify "BLOB ENTRY lsq $LSQ: λ₁=$L1"
            echo "Blob entry at lsq $LSQ"
            ZONE_STATE="blob"
        fi
    elif ! echo "$line" | grep -q "<--"; then
        if echo "$line" | grep -qE '^ [0-9]'; then
            ZONE_STATE="blob"
        fi
    fi
done
