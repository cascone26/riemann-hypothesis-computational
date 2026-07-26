#!/usr/bin/env python3
"""
Phase 2: Attempt to extend toward ½log5 ≈ 0.8047.
Use width 0.04 for larger exploratory windows.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import json
from datetime import datetime
from sos_affine_compose import build_affine_entry_enclosures, sos_certify_from_intervals

results_phase2 = []
a_extend_start = 0.55
half_log5 = 0.8047
width = 0.04

# Create windows
windows = []
a = a_extend_start
while a < half_log5:
    a_next = min(a + width, half_log5)
    windows.append((a, a_next))
    a = a_next

print("="*100)
print(f"PHASE 2: Extend [0.55, {half_log5:.4f}] with width {width}")
print("="*100)

for i, (a0, a1) in enumerate(windows):
    w = a1 - a0
    print(f"[{i+1}/{len(windows)}] [{a0:.4f}, {a1:.4f}]...", end=" ", flush=True)

    try:
        lo_mat, hi_mat = build_affine_entry_enclosures(a0, a1)
        widths = hi_mat - lo_mat
        status, cert = sos_certify_from_intervals(a0, a1, lo_mat, hi_mat)

        if cert is None:
            print(f"FAIL ({status})")
            results_phase2.append({
                "a0": float(a0), "a1": float(a1), "width": float(w),
                "status": str(status), "gram_min_eig": None, "passed": False
            })
        else:
            A, B, C, P1 = cert
            Gram = np.block([[A, B], [B.T, C]])
            gmin = np.linalg.eigvalsh(Gram).min()
            pmin = np.linalg.eigvalsh(P1).min()
            passed = gmin >= -1e-6 and pmin >= -1e-6

            print(f"{'PASS' if passed else 'FAIL'} (gram={gmin:.2e})")
            results_phase2.append({
                "a0": float(a0), "a1": float(a1), "width": float(w),
                "status": str(status), "gram_min_eig": float(gmin),
                "entry_width_max": float(widths.max()), "passed": bool(passed)
            })
    except Exception as e:
        print(f"ERROR: {e}")
        results_phase2.append({
            "a0": float(a0), "a1": float(a1), "width": float(w),
            "error": str(e), "passed": False
        })

print("\nPHASE 2 SUMMARY:")
passed = len([r for r in results_phase2 if r.get('passed', False)])
failed = len([r for r in results_phase2 if not r.get('passed', False)])
print(f"  Total: {len(results_phase2)}, Passed: {passed}, Failed: {failed}")

if failed:
    print("  Failed windows:")
    for r in results_phase2:
        if not r.get('passed', False):
            print(f"    [{r['a0']:.4f}, {r['a1']:.4f}]")

# Save results
report = {
    "phase": 2,
    "interval": [a_extend_start, half_log5],
    "window_width": width,
    "windows": len(results_phase2),
    "passed": passed,
    "failed": failed,
    "results": results_phase2,
    "timestamp": datetime.now().isoformat() + "Z"
}

os.makedirs("../../results/track_continuum", exist_ok=True)
with open("../../results/track_continuum/phase2_results.json", "w") as f:
    json.dump(report, f, indent=2)

print(f"\nResults saved to: results/track_continuum/phase2_results.json")
print("="*100)
