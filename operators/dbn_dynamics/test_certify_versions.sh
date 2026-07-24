#!/bin/bash
# Test harness: run all optimization versions and compare

cd /Users/scones/projects/riemann

PYTHON=/Users/scones/projects/riemann/flint_venv/bin/python3

echo "=== CERTIFICATION SCRIPT OPTIMIZATION TEST ==="
echo "Running all versions and comparing runtimes..."
echo

# Helper function to extract timing and result
extract_result() {
    local file=$1
    local name=$2
    if [ ! -f "$file" ]; then
        echo "$name: NOT FOUND"
        return
    fi
    if ! grep -q "CERTIFIED\|not certified" "$file"; then
        echo "$name: INCOMPLETE (still running or error)"
        return
    fi
    local runtime=$(grep "runtime=" "$file" | head -1 | sed 's/.*runtime=\([0-9.]*\)s.*/\1/')
    if [ -z "$runtime" ]; then
        runtime=$(grep "elapsed:" "$file" | head -1 | sed 's/.*elapsed: \([0-9.]*\)s.*/\1/')
    fi
    local certified=$(grep -q "*** CERTIFIED" "$file" && echo "PASS" || echo "FAIL")
    echo "$name: ${runtime:----}s $certified"
    grep "pivot" "$file" | head -4
    echo
}

# Wait for all files to have results
echo "Waiting for all tests to complete..."
until [ -f /tmp/orig_test.txt ] && [ -f /tmp/opt_test.txt ] && [ -f /tmp/opt2_test.txt ] && [ -f /tmp/opt3_test.txt ] && \
      grep -q "CERTIFIED\|not certified" /tmp/orig_test.txt 2>/dev/null && \
      grep -q "CERTIFIED\|not certified" /tmp/opt_test.txt 2>/dev/null && \
      grep -q "CERTIFIED\|not certified" /tmp/opt2_test.txt 2>/dev/null && \
      grep -q "CERTIFIED\|not certified" /tmp/opt3_test.txt 2>/dev/null; do
    sleep 5
    echo -n "."
done
echo " Done!"
echo

echo "=== RESULTS ==="
echo
extract_result "/tmp/orig_test.txt" "ORIGINAL v2"
extract_result "/tmp/opt_test.txt"  "OPT v1 (250-bit precision)"
extract_result "/tmp/opt2_test.txt" "OPT v2 (context switch)"
extract_result "/tmp/opt3_test.txt" "OPT v3 (Lerch constant split)"

echo "=== FULL OUTPUTS ==="
echo
echo "--- ORIGINAL ---"
cat /tmp/orig_test.txt
echo
echo "--- OPT v1 ---"
cat /tmp/opt_test.txt
echo
echo "--- OPT v2 ---"
cat /tmp/opt2_test.txt
echo
echo "--- OPT v3 ---"
cat /tmp/opt3_test.txt
