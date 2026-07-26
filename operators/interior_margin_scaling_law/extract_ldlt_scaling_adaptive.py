#!/usr/bin/env python3
"""Extract interior margin (min pivot) vs conductor c scaling law with adaptive precision.

Runs arb_ldlt_certify.py for c in a specified range with adaptive precision
that increases with c to maintain certification accuracy.
"""

import sys
import os
import subprocess
import json
import tempfile
import argparse
from pathlib import Path
import math

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


def certified_inertia_with_pivots(A, DIM, heartbeat=None):
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
            pct = 100.0 * (i + 1) / DIM
            print(f"    {pct:5.1f}% pivot {i + 1}/{DIM}  elapsed={el:.0f}s", flush=True)
    return n_pos, n_neg, None, d


def extract_min_pivot(pivots):
    """Extract minimum certified pivot value (absolute distance from 0)."""
    min_abs_pivot = None
    for p in pivots:
        if p is None:
            continue
        try:
            # For arb balls, get the lower bound of absolute value
            if hasattr(p, 'a'):
                left = float(p.a)
                right = float(p.b)
            else:
                left = right = float(p)

            # Distance to 0 is minimum of |left| and |right|
            abs_dist_to_0 = min(abs(left), abs(right))
        except:
            continue

        if min_abs_pivot is None or abs_dist_to_0 < min_abs_pivot:
            min_abs_pivot = abs_dist_to_0

    return min_abs_pivot


def get_precision_for_c(c, base_prec=600):
    """Return precision bits needed for conductor c.

    Empirically, precision needs to grow roughly logarithmically with c,
    to keep the minimum pivot certified accurately.
    """
    # Based on testing: c=2-25 work fine with 600-800
    # c=30+ need 1200+, c=50+ need 2000+, etc.
    if c <= 20:
        return base_prec
    elif c <= 30:
        return int(base_prec * 1.5)
    elif c <= 50:
        return int(base_prec * 2.5)
    elif c <= 100:
        return int(base_prec * 4)
    else:
        return int(base_prec * 5)


def run_one_conductor(c, N, prec, max_retries=2):
    """Run certifier for one conductor value. Returns (c, λ_min, n_pos, n_neg) or None."""
    print(f"[RUN] c={c}, N={N}, prec={prec}", flush=True)

    for attempt in range(max_retries):
        try:
            if attempt > 0:
                prec_attempt = int(prec * (1.5 ** attempt))
                print(f"  Retry {attempt} with prec={prec_attempt}", flush=True)
            else:
                prec_attempt = prec

            # Build the matrix
            A, DIM = certify_module.build_arb_tau(c, N, prec_attempt)

            # Run LDL^T with heartbeat every 20 pivots for progress
            heartbeat = max(20, DIM // 10)
            n_pos, n_neg, undet, pivots = certified_inertia_with_pivots(A, DIM, heartbeat=heartbeat)

            if undet is not None:
                if attempt < max_retries - 1:
                    continue  # Retry with higher precision
                else:
                    print(f"  UNDETERMINED at pivot {undet} even after retries", flush=True)
                    return None

            # Extract λ_min
            lambda_min = extract_min_pivot(pivots)
            if lambda_min is None:
                print(f"  ERROR: Could not extract λ_min", flush=True)
                return None

            print(f"  ✓ λ_min = {lambda_min:.6e}  n_pos={n_pos}, n_neg={n_neg}", flush=True)
            return c, lambda_min, n_pos, n_neg, DIM

        except Exception as e:
            if attempt < max_retries - 1:
                print(f"  Error on attempt {attempt+1}: {e}, retrying...", flush=True)
            else:
                print(f"  ERROR: {e}", flush=True)
                return None

    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-start", type=int, default=2, help="Start conductor")
    parser.add_argument("--c-stop", type=int, default=100, help="Stop conductor")
    parser.add_argument("--c-skip-above", type=int, default=None, help="Skip c > this value")
    parser.add_argument("--N", type=int, default=100, help="Truncation degree N")
    parser.add_argument("--base-prec", type=int, default=600, help="Base Arb precision in bits")
    parser.add_argument("--output", type=str, default="scaling_law_adaptive.json")
    args = parser.parse_args()

    results = []
    failed = []

    skip_above = args.c_skip_above if args.c_skip_above else args.c_stop

    for c in range(args.c_start, min(args.c_stop, skip_above) + 1):
        prec = get_precision_for_c(c, args.base_prec)
        try:
            result = run_one_conductor(c, args.N, prec)
            if result is not None:
                c_val, lambda_min, n_pos, n_neg, dim = result
                results.append({
                    "c": c_val,
                    "lambda_min": lambda_min,
                    "log10_lambda_min": float(mp.log(lambda_min, 10)) if lambda_min > 0 else None,
                    "log_lambda_min": float(mp.log(lambda_min)) if lambda_min > 0 else None,
                    "n_pos": n_pos,
                    "n_neg": n_neg,
                    "dimension": dim,
                    "N": args.N,
                    "prec_used": prec
                })
            else:
                failed.append(c)
        except KeyboardInterrupt:
            print(f"\nInterrupted at c={c}")
            break
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
                "c_skip_above": skip_above,
                "N": args.N,
                "base_prec": args.base_prec,
                "failed": failed,
                "count": len(results)
            },
            "data": results
        }, f, indent=2)

    print(f"\n✓ Saved {len(results)} results to {output_path}")
    if failed:
        print(f"  Failed at c: {failed}")

    # Print summary
    if results:
        print("\n=== SCALING LAW DATA ===")
        print("c    | λ_min        | log10(λ_min) | log(λ_min)")
        print("-----|--------------|--------------|----------")
        for r in results:
            lmin_str = f"{r['lambda_min']:.6e}" if r['lambda_min'] else "None"
            log10_str = f"{r['log10_lambda_min']:.2f}" if r['log10_lambda_min'] else "None"
            log_str = f"{r['log_lambda_min']:.2f}" if r['log_lambda_min'] else "None"
            print(f"{r['c']:4d} | {lmin_str:12s} | {log10_str:>12s} | {log_str:>8s}")


if __name__ == "__main__":
    main()
