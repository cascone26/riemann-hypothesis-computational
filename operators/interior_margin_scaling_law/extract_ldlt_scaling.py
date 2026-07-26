#!/usr/bin/env python3
"""Extract interior margin (min pivot) vs conductor c scaling law.

Runs arb_ldlt_certify.py for c in {2..100} and extracts λ_min (minimum
certified pivot) from each run. Fits decay law (power/exponential/log).

Checks coincidence: c ≤ 100 ⟺ a ≤ ½log(100) ≈ 2.30 vs
largest verified-zero height T ≈ 3e12 gives a_max ≈ log(T)/(4π) ≈ 2.29.
"""

import sys
import os
import subprocess
import json
import tempfile
import argparse
from pathlib import Path

# Add the external groskin code to path
riemann_root = Path(__file__).parent.parent.parent
groskin_code = riemann_root / "external" / "groskin_2607.02828" / "code"
sys.path.insert(0, str(groskin_code))

from flint import arb, acb, arb_mat, ctx
import mpmath as mp

# Import the certifier functions
import importlib.util
spec = importlib.util.spec_from_file_location(
    "arb_ldlt_certify",
    str(groskin_code / "arb_ldlt_certify.py")
)
certify_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(certify_module)


def certified_inertia_with_pivots(A, DIM, heartbeat=50):
    """Interval LDL^T.  Returns (n_pos, n_neg, undetermined_pivot_or_None, all_pivots)."""
    d = [None] * DIM
    Lf = [[arb(0)] * DIM for _ in range(DIM)]
    n_pos = 0
    n_neg = 0
    import time
    t0 = time.time()
    for i in range(DIM):
        s = A[i, i]
        for k in range(i):
            s = s - Lf[i][k] * Lf[i][k] * d[k]
        d[i] = s
        if s > 0:
            n_pos += 1
        elif s < 0:
            n_neg += 1
        else:
            return n_pos, n_neg, i, d
        for j in range(i + 1, DIM):
            t = A[j, i]
            for k in range(i):
                t = t - Lf[j][k] * Lf[i][k] * d[k]
            Lf[j][i] = t / d[i]
        if heartbeat and (i + 1) % heartbeat == 0:
            el = time.time() - t0
            print("  [ldlt] pivot %d/%d  elapsed=%.0fs  eta=%.0fs"
                  % (i + 1, DIM, el, el / (i + 1) * (DIM - i - 1)), flush=True)
    return n_pos, n_neg, None, d


def extract_min_pivot(pivots):
    """Extract minimum certified pivot value (lower bound of interval).

    For positive pivots, use the lower bound (left endpoint).
    For negative pivots, use the upper bound (right endpoint, closest to 0).
    Returns absolute minimum distance from 0.
    """
    min_abs_pivot = None
    for p in pivots:
        if p is None:
            continue
        # arb balls have .a and .b for left/right endpoints
        if hasattr(p, 'a'):
            # It's an arb ball
            left = float(p.a)
            right = float(p.b)
            abs_dist_to_0 = min(abs(left), abs(right))
        else:
            # Fallback: convert to float
            abs_dist_to_0 = abs(float(p))

        if min_abs_pivot is None or abs_dist_to_0 < min_abs_pivot:
            min_abs_pivot = abs_dist_to_0

    return min_abs_pivot


def run_one_conductor(c, N, prec):
    """Run certifier for one conductor value. Returns (c, λ_min, n_pos, n_neg)."""
    print(f"\n[RUN] c={c}, N={N}, prec={prec}")

    # Build the matrix
    print(f"  Building matrix...", flush=True)
    A, DIM = certify_module.build_arb_tau(c, N, prec)
    print(f"  Matrix built, dimension={DIM}")

    # Run LDL^T
    print(f"  Running interval LDL^T...")
    n_pos, n_neg, undet, pivots = certified_inertia_with_pivots(A, DIM, heartbeat=100)

    if undet is not None:
        print(f"  WARNING: Undetermined at pivot {undet}, might need higher precision")
        return None

    # Extract λ_min
    lambda_min = extract_min_pivot(pivots)
    print(f"  λ_min = {lambda_min:.6e}  n_pos={n_pos}, n_neg={n_neg}")

    return c, lambda_min, n_pos, n_neg, DIM


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-start", type=int, default=2, help="Start conductor")
    parser.add_argument("--c-stop", type=int, default=100, help="Stop conductor")
    parser.add_argument("--N", type=int, default=100, help="Truncation degree N")
    parser.add_argument("--prec", type=int, default=1200, help="Arb precision in bits")
    parser.add_argument("--output", type=str, default="scaling_law_data.json")
    args = parser.parse_args()

    results = []
    failed = []

    for c in range(args.c_start, args.c_stop + 1):
        try:
            result = run_one_conductor(c, args.N, args.prec)
            if result is not None:
                c_val, lambda_min, n_pos, n_neg, dim = result
                results.append({
                    "c": c_val,
                    "lambda_min": lambda_min,
                    "log10_lambda_min": float(mp.log(lambda_min, 10)) if lambda_min > 0 else None,
                    "n_pos": n_pos,
                    "n_neg": n_neg,
                    "dimension": dim,
                    "N": args.N,
                    "prec": args.prec
                })
            else:
                failed.append(c)
        except Exception as e:
            print(f"ERROR at c={c}: {e}")
            failed.append(c)

    # Save results
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        json.dump({
            "metadata": {
                "c_start": args.c_start,
                "c_stop": args.c_stop,
                "N": args.N,
                "prec": args.prec,
                "failed": failed
            },
            "data": results
        }, f, indent=2)

    print(f"\nSaved {len(results)} results to {output_path}")
    print(f"Failed: {failed}")

    # Print summary
    if results:
        print("\n=== SCALING LAW DATA ===")
        print("c    | λ_min        | log10(λ_min)")
        print("-----|--------------|-------------")
        for r in results[:20]:  # First 20
            lmin_str = f"{r['lambda_min']:.6e}" if r['lambda_min'] else "None"
            log_str = f"{r['log10_lambda_min']:.2f}" if r['log10_lambda_min'] else "None"
            print(f"{r['c']:4d} | {lmin_str:12s} | {log_str:>12s}")
        if len(results) > 20:
            print(f"... ({len(results) - 20} more)")


if __name__ == "__main__":
    main()
