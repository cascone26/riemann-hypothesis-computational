#!/bin/bash
# Quick launcher for c=100 analysis
# Run from riemann project root: ./operators/RUN_C100_ANALYSIS.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
SCRIPT="$SCRIPT_DIR/analyze_c100.py"

echo "======================================================================="
echo "CvS c=100 Analysis Launcher"
echo "======================================================================="
echo ""
echo "Project root: $PROJECT_ROOT"
echo "Script: $SCRIPT"
echo ""

# Check script exists
if [[ ! -f "$SCRIPT" ]]; then
    echo "ERROR: analyze_c100.py not found at $SCRIPT"
    exit 1
fi

# Check mpmath
if ! python3 -c "import mpmath" 2>/dev/null; then
    echo "WARNING: mpmath not installed. Installing..."
    pip install mpmath
fi

# Check SSH access to Olympus
echo "Checking SSH access to olympus..."
if ! ssh olympus echo "OK" 2>/dev/null; then
    echo "ERROR: Cannot reach olympus via SSH"
    echo ""
    echo "Troubleshooting:"
    echo "  1. Check Tailscale: ping olympus"
    echo "  2. Check SSH key: ssh -v olympus"
    echo "  3. Ensure SSH config has olympus alias (or use full IP)"
    exit 1
fi
echo "✓ SSH access OK"
echo ""

# Check remote result file
echo "Checking if remote result file exists..."
if ssh olympus test -f /home/cobo/riemann/results/cvs_highprec_c100.json; then
    SIZE=$(ssh olympus wc -c < /home/cobo/riemann/results/cvs_highprec_c100.json | tr -d ' ')
    echo "✓ Result file exists ($(numfmt --to=iec $SIZE 2>/dev/null || echo $SIZE bytes))"
else
    echo "WARNING: Result file not found on olympus yet"
    echo "  Path: /home/cobo/riemann/results/cvs_highprec_c100.json"
    echo "  The computation may still be running."
    echo ""
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi
echo ""

# Run analysis
echo "======================================================================="
echo "Running analysis..."
echo "======================================================================="
cd "$PROJECT_ROOT"
python3 "$SCRIPT"

echo ""
echo "======================================================================="
echo "Done! Check results in:"
echo "  $PROJECT_ROOT/results/SESSION_FINDINGS_2026_07_11.md"
echo "======================================================================="
