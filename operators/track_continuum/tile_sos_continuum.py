#!/usr/bin/env python3
"""
Tile the interval [0.35, 0.55] (and optionally extend toward 0.8047) with adjacent
SOS certification windows. Report which windows certify Q(a)⪰0 and which fail.

REPORT TIER: Float64 SDP solver output, not re-verified with rational rounding or
extended-precision. See sos_affine_compose.py header for Handle-tier notes.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'operators', 'track_continuum'))

import numpy as np
import json
from datetime import datetime

# Import the SOS certifier functions
from sos_affine_compose import build_affine_entry_enclosures, sos_certify_from_intervals

def tile_and_certify(windows, verbose=True):
    """
    Run SOS certification on a list of (a0, a1) windows.
    Returns a list of dicts: {a0, a1, width, status, gram_min_eig, passed}
    """
    results = []

    for i, (a0, a1) in enumerate(windows):
        width = a1 - a0
        if verbose:
            print(f"\n[{i+1}/{len(windows)}] Window [{a0:.4f}, {a1:.4f}] (width {width:.4f})")

        # Build affine enclosures
        lo_mat, hi_mat = build_affine_entry_enclosures(a0, a1)
        widths = hi_mat - lo_mat

        # SOS certify
        status, cert = sos_certify_from_intervals(a0, a1, lo_mat, hi_mat)

        if cert is None:
            if verbose:
                print(f"  Status: {status} (FAIL - no feasible certificate)")
            results.append({
                "a0": float(a0),
                "a1": float(a1),
                "width": float(width),
                "status": str(status),
                "gram_min_eig": None,
                "entry_width_max": float(widths.max()),
                "passed": False
            })
        else:
            import numpy as np
            A, B, C, P1 = cert
            Gram = np.block([[A, B], [B.T, C]])
            gmin = np.linalg.eigvalsh(Gram).min()
            pmin = np.linalg.eigvalsh(P1).min()

            passed = gmin >= -1e-6 and pmin >= -1e-6
            if verbose:
                print(f"  Status: {status}")
                print(f"    Gram min eig: {gmin:.3e}, P1 min eig: {pmin:.3e}")
                print(f"    Entry width max: {widths.max():.3e}")
                print(f"  Result: {'PASS' if passed else 'FAIL'}")

            results.append({
                "a0": float(a0),
                "a1": float(a1),
                "width": float(width),
                "status": str(status),
                "gram_min_eig": float(gmin),
                "entry_width_max": float(widths.max()),
                "passed": passed
            })

    return results

if __name__ == "__main__":
    print("="*100)
    print("CONTINUUM TILING: SOS certification across windowed intervals")
    print("Target: p=2-only interval [½log2, ½log3] ≈ [0.3466, 0.5493]")
    print("Task: tile [0.35, 0.55] with windows of width 0.03")
    print("REPORT TIER: Float64 SDP (not Handle-tier; see sos_affine_compose.py)")
    print("="*100)

    # Phase 1: Main tiling [0.35, 0.55] with width 0.03
    print("\n[PHASE 1] Tiling [0.35, 0.55] with width 0.03")
    a_start, a_end = 0.35, 0.55
    width = 0.03
    n_windows = int(np.ceil((a_end - a_start) / width))
    windows_1 = [(a_start + i*width, a_start + (i+1)*width) for i in range(n_windows)]

    # Adjust last window to end at a_end
    if windows_1[-1][1] > a_end:
        windows_1[-1] = (windows_1[-1][0], a_end)

    print(f"Creating {len(windows_1)} windows:")
    for a0, a1 in windows_1:
        print(f"  [{a0:.4f}, {a1:.4f}]")

    results_1 = tile_and_certify(windows_1, verbose=True)

    # Summarize Phase 1
    passed = [r for r in results_1 if r['passed']]
    failed = [r for r in results_1 if not r['passed']]

    print(f"\n[PHASE 1 SUMMARY]")
    print(f"  Total windows: {len(results_1)}")
    print(f"  Passed: {len(passed)}")
    print(f"  Failed: {len(failed)}")

    if failed:
        print(f"  Failed windows:")
        for r in failed:
            print(f"    [{r['a0']:.4f}, {r['a1']:.4f}]: {r['status']}, gram_min={r['gram_min_eig']}")

    # Phase 2: Attempt to extend toward ½log5 ≈ 0.8047
    print("\n[PHASE 2] Attempt to extend toward ½log5 ≈ 0.8047")
    half_log5 = 0.8047
    a_extend_start = a_end  # 0.55

    # Use larger windows (0.04-0.05) for this exploratory phase
    width_2 = 0.04
    n_windows_2 = int(np.ceil((half_log5 - a_extend_start) / width_2))
    windows_2 = [(a_extend_start + i*width_2, a_extend_start + (i+1)*width_2) for i in range(n_windows_2)]

    # Adjust last window
    if windows_2[-1][1] > half_log5:
        windows_2[-1] = (windows_2[-1][0], half_log5)

    print(f"Creating {len(windows_2)} extension windows:")
    for a0, a1 in windows_2:
        print(f"  [{a0:.4f}, {a1:.4f}]")

    results_2 = tile_and_certify(windows_2, verbose=True)

    # Summarize Phase 2
    passed_2 = [r for r in results_2 if r['passed']]
    failed_2 = [r for r in results_2 if not r['passed']]

    print(f"\n[PHASE 2 SUMMARY]")
    print(f"  Total windows: {len(results_2)}")
    print(f"  Passed: {len(passed_2)}")
    print(f"  Failed: {len(failed_2)}")

    if failed_2:
        print(f"  Failed windows:")
        for r in failed_2:
            print(f"    [{r['a0']:.4f}, {r['a1']:.4f}]: {r['status']}, gram_min={r['gram_min_eig']}")

    # Save consolidated report
    all_results = results_1 + results_2
    report = {
        "tool": "tile_sos_continuum.py",
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "tier": "REPORT (float64 SDP, not extended-precision verified)",
        "theorem": "THEOREM_ASSEMBLY piece 6: continuum tiling certification",
        "intervals": {
            "phase1": {
                "range": [a_start, a_end],
                "window_width": width,
                "windows": len(windows_1),
                "passed": len(passed),
                "failed": len(failed),
                "results": results_1
            },
            "phase2": {
                "range": [a_extend_start, half_log5],
                "window_width": width_2,
                "windows": len(windows_2),
                "passed": len(passed_2),
                "failed": len(failed_2),
                "results": results_2
            }
        },
        "summary": {
            "total_windows": len(all_results),
            "total_passed": len(passed) + len(passed_2),
            "total_failed": len(failed) + len(failed_2),
            "coverage": {
                "p2_only": f"[0.35, 0.55] fully tiled" if len(failed) == 0 else f"gaps detected",
                "extended": f"[0.55, {half_log5:.4f}] {'fully tiled' if len(failed_2) == 0 else 'partial coverage'}"
            }
        }
    }

    os.makedirs("../../results/track_continuum", exist_ok=True)
    report_file = "../../results/track_continuum/tile_sos_continuum_report.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)

    print(f"\n[FINAL REPORT]")
    print(f"  Total windows certified: {len(passed) + len(passed_2)} / {len(all_results)}")
    print(f"  Report saved: {report_file}")
    print("="*100)
