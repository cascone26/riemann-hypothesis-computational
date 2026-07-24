#!/usr/bin/env python3
"""Wrapper for ODD sector certification. Calls arb_certify.py --parity odd.
This allows verify_gate.py's "bound" type (which can't pass arguments) to
re-verify the odd-sector cert independently."""
import subprocess
import sys
result = subprocess.run([
    sys.executable,
    '/Users/scones/projects/riemann/operators/dbn_dynamics/arb_certify.py',
    '--parity', 'odd'
], capture_output=False, text=True)
sys.exit(result.returncode)
