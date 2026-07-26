#!/usr/bin/env python3
"""
Phase 1: Tile [0.35, 0.55] with width 0.03.
Simpler version focused on getting results.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import json
from datetime import datetime
from sos_affine_compose import build_affine_entry_enclosures, sos_certify_from_intervals

results_phase1 = []
windows = [
    (0.35, 0.38),
    (0.38, 0.41),
    (0.41, 0.44),
    (0.44, 0.47),
    (0.47, 0.50),
    (0.50, 0.53),
    (0.53, 0.55),
]

print("="*100)
print("PHASE 1: Tiling [0.35, 0.55] with width 0.03")
print("="*100)

for i, (a0, a1) in enumerate(windows):
    width = a1 - a0
    print(f"[{i+1}/{len(windows)}] [{a0:.4f}, {a1:.4f}]...", end=" ", flush=True)

    try:
        lo_mat, hi_mat = build_affine_entry_enclosures(a0, a1)
        widths = hi_mat - lo_mat
        status, cert = sos_certify_from_intervals(a0, a1, lo_mat, hi_mat)

        if cert is None:
            print(f"FAIL ({status})")
            results_phase1.append({
                "a0": float(a0), "a1": float(a1), "width": float(width),
                "status": str(status), "gram_min_eig": None, "passed": False
            })
        else:
            A, B, C, P1 = cert
            Gram = np.block([[A, B], [B.T, C]])
            gmin = np.linalg.eigvalsh(Gram).min()
            pmin = np.linalg.eigvalsh(P1).min()
            passed = gmin >= -1e-6 and pmin >= -1e-6

            print(f"{'PASS' if passed else 'FAIL'} (gram={gmin:.2e})")
            results_phase1.append({
                "a0": float(a0), "a1": float(a1), "width": float(width),
                "status": str(status), "gram_min_eig": float(gmin),
                "entry_width_max": float(widths.max()), "passed": passed
            })
    except Exception as e:
        print(f"ERROR: {e}")
        results_phase1.append({
            "a0": float(a0), "a1": float(a1), "width": float(width),
            "error": str(e), "passed": False
        })

print("\nPHASE 1 SUMMARY:")
passed = len([r for r in results_phase1 if r['passed']])
failed = len([r for r in results_phase1 if not r['passed']])
print(f"  Total: {len(results_phase1)}, Passed: {passed}, Failed: {failed}")

if failed:
    print("  Failed windows:")
    for r in results_phase1:
        if not r['passed']:
            print(f"    [{r['a0']:.4f}, {r['a1']:.4f}]")

# Save results
report = {
    "phase": 1,
    "interval": [0.35, 0.55],
    "window_width": 0.03,
    "windows": len(results_phase1),
    "passed": passed,
    "failed": failed,
    "results": results_phase1,
    "timestamp": datetime.utcnow().isoformat() + "Z"
}

os.makedirs("../../results/track_continuum", exist_ok=True)
with open("../../results/track_continuum/phase1_results.json", "w") as f:
    json.dump(report, f, indent=2)

print(f"\nResults saved to: results/track_continuum/phase1_results.json")
print("="*100)
