#!/usr/bin/env python3
"""
Quick test: run SOS certification on just 3 windows to validate the tiling approach.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'operators', 'track_continuum'))

import numpy as np
import json
from datetime import datetime

# Import the SOS certifier functions
from sos_affine_compose import build_affine_entry_enclosures, sos_certify_from_intervals

def test_windows():
    """Test on 3 quick windows."""
    windows = [
        (0.35, 0.38),
        (0.45, 0.48),
        (0.55, 0.58)
    ]

    results = []
    for i, (a0, a1) in enumerate(windows):
        width = a1 - a0
        print(f"Window {i+1}: [{a0}, {a1}]", flush=True)
        sys.stdout.flush()

        # Build affine enclosures
        try:
            lo_mat, hi_mat = build_affine_entry_enclosures(a0, a1)
            widths = hi_mat - lo_mat
            print(f"  Affine enclosures built", flush=True)
            sys.stdout.flush()
        except Exception as e:
            print(f"  ERROR in affine: {e}", flush=True)
            sys.stdout.flush()
            continue

        # SOS certify
        try:
            status, cert = sos_certify_from_intervals(a0, a1, lo_mat, hi_mat)
            print(f"  SOS status: {status}", flush=True)
            sys.stdout.flush()

            if cert is None:
                results.append({"a0": a0, "a1": a1, "passed": False, "status": status})
                print(f"  Result: FAIL", flush=True)
            else:
                A, B, C, P1 = cert
                Gram = np.block([[A, B], [B.T, C]])
                gmin = np.linalg.eigvalsh(Gram).min()
                passed = gmin >= -1e-6
                results.append({"a0": a0, "a1": a1, "passed": passed, "gram_min": float(gmin)})
                print(f"  Result: {'PASS' if passed else 'FAIL'} (gram_min={gmin:.3e})", flush=True)
            sys.stdout.flush()
        except Exception as e:
            print(f"  ERROR in SOS: {e}", flush=True)
            import traceback
            traceback.print_exc()
            sys.stdout.flush()

    passed = len([r for r in results if r['passed']])
    print(f"\nSummary: {passed}/{len(results)} windows passed")
    return results

if __name__ == "__main__":
    print("Testing tile_sos_continuum with 3 windows", flush=True)
    sys.stdout.flush()
    test_windows()
