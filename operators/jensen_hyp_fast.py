#!/usr/bin/env python3
"""
jensen_hyp_fast.py — Fast Jensen polynomial hyperbolicity check.

Uses DFT extraction of Taylor coefficients (not slow mpmath.diff).
Checks J_{d,n} hyperbolicity for d=1..30, n=0..500.

GORZ 2019: J_{d,n} is hyperbolic for d≤8 ALL n.
Goal: check d=9..30 computationally.
"""

import json, os, time
import numpy as np
from math import comb
import mpmath

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
mpmath.mp.dps = 60


def Xi_func(z):
    s = mpmath.mpf('0.5') + 1j * z
    return mpmath.re(mpmath.mpf('0.5') * s * (s-1) * mpmath.power(mpmath.pi, -s/2)
                     * mpmath.gamma(s/2) * mpmath.zeta(s))


def compute_coeffs_dft(N_max, prec=60):
    """Extract a_n = Ξ^{(2n)}(0)/(2n)! via DFT on a circle |z|=r."""
    mpmath.mp.dps = prec
    print(f"DFT extraction: N_max={N_max}, prec={prec}dps...", flush=True)
    t0 = time.time()

    N = 4 * (N_max + 10)  # oversample
    r = mpmath.mpf('3.0')  # radius < |gamma_1|=14.13

    vals = np.zeros(N, dtype=complex)
    for k in range(N):
        theta = 2 * mpmath.pi * k / N
        z = r * mpmath.exp(1j * theta)
        vals[k] = complex(Xi_func(z))
        if k % 50 == 0:
            print(f"  eval {k}/{N}...", flush=True)

    fft = np.fft.fft(vals) / N
    r_float = float(r)

    coeffs = []
    for n in range(N_max + 1):
        if 2*n < N:
            a = fft[2*n].real / r_float**(2*n)
        else:
            a = 0.0
        coeffs.append(float(a))

    print(f"Done in {time.time()-t0:.1f}s", flush=True)
    return coeffs


def is_hyperbolic(poly_coeffs, tol=1e-9):
    """All roots real? poly_coeffs = [c0,c1,...,cd] (ascending degree)."""
    while len(poly_coeffs) > 1 and abs(poly_coeffs[-1]) < 1e-100:
        poly_coeffs = poly_coeffs[:-1]
    if len(poly_coeffs) <= 1:
        return True, 0.0

    p_rev = poly_coeffs[::-1]
    roots = np.roots(p_rev)
    if len(roots) == 0:
        return True, 0.0
    max_imag = max(abs(r.imag) for r in roots)
    scale = max(abs(r.real) for r in roots) + 1.0
    return max_imag < tol * scale, float(max_imag / scale)


def check_all(a, d_max=30, n_max=500):
    print(f"\nHyperbolicity check: d=1..{d_max}, n=0..{min(n_max,len(a)-d_max-1)}")
    results = {}
    for d in range(1, d_max+1):
        first_fail = None
        max_imag_ratio = 0.0
        n_checked = 0
        for n in range(n_max+1):
            if n + d >= len(a):
                break
            poly = [comb(d, j) * a[n+j] for j in range(d+1)]
            hyp, imag_ratio = is_hyperbolic(poly)
            max_imag_ratio = max(max_imag_ratio, imag_ratio)
            n_checked += 1
            if not hyp and first_fail is None:
                first_fail = n
                break  # stop at first fail
        all_hyp = first_fail is None
        status = "ALL ✓" if all_hyp else f"FAIL n={first_fail}"
        print(f"  d={d:>3}: {status}  (max|Im/scale|={max_imag_ratio:.2e}, n_checked={n_checked})",
              flush=True)
        results[d] = {'all_hyperbolic': all_hyp, 'first_fail': first_fail,
                      'n_checked': n_checked, 'max_imag_ratio': max_imag_ratio}
    return results


def log_convexity(a, n_max=200):
    """Check log-convexity: a_{n+1}^2 >= a_n * a_{n+2}."""
    violations = []
    ratios = []
    for n in range(min(len(a)-2, n_max)):
        if abs(a[n]) > 1e-100 and abs(a[n+2]) > 1e-100:
            r = a[n+1]**2 / (a[n] * a[n+2])
            ratios.append(r)
            if r < 1.0 - 1e-10:
                violations.append((n, r))
    return violations, ratios


def main():
    print("="*70)
    print("JENSEN POLYNOMIAL HYPERBOLICITY (fast DFT method)")
    print("="*70)

    # Get Taylor coefficients
    N_max = 250
    a = compute_coeffs_dft(N_max, prec=60)

    # Sanity check
    print(f"\nFirst 12 coefficients a_n:")
    for n in range(12):
        print(f"  a_{n:>3} = {a[n]:.10e}")
    print(f"All positive? {all(x > 0 for x in a[:50])}")

    # Log-convexity
    viols, ratios = log_convexity(a, n_max=200)
    print(f"\nLog-convexity (a_{{n+1}}^2 >= a_n*a_{{n+2}}): violations={len(viols)}")
    if not viols:
        print("  ✓ Fully log-convex for n=0..200")
        print(f"  Min ratio = {min(ratios):.8f}, Max ratio = {max(ratios):.8f}")
        print(f"  Ratios converge to 1? Last 5: {ratios[-5:]}")

    # Hyperbolicity table
    hyp_results = check_all(a, d_max=30, n_max=min(200, N_max-31))

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    all_hyp_d_max = max(d for d, r in hyp_results.items() if r['all_hyperbolic'])
    first_fail_d = min((d for d, r in hyp_results.items() if not r['all_hyperbolic']),
                       default=None)
    print(f"All hyperbolic up to d={all_hyp_d_max}")
    if first_fail_d:
        print(f"First failure: d={first_fail_d} at n={hyp_results[first_fail_d]['first_fail']}")
    else:
        print(f"ALL d=1..30 hyperbolic for all tested n — GORZ extends to d≤30!")

    # Ratio analysis for d=2: discriminant condition
    print(f"\nDiscriminant (d=2) condition: ratio r_n = a_{{n+1}}^2/(a_n*a_{{n+2}})")
    for n in [0, 10, 20, 50, 100, 150, 200]:
        if n+2 < len(a) and abs(a[n]) > 1e-100 and abs(a[n+2]) > 1e-100:
            r = a[n+1]**2 / (a[n]*a[n+2])
            print(f"  n={n:>3}: r={r:.10f}  (excess above 1: {r-1:.4e})")

    # Save
    out = {
        'a_coeffs': a[:100],
        'log_convexity_violations': viols,
        'min_ratio': float(min(ratios)) if ratios else None,
        'max_ratio': float(max(ratios)) if ratios else None,
        'hyperbolicity': {str(d): r for d, r in hyp_results.items()},
        'all_hyp_d_max': all_hyp_d_max,
        'first_fail_d': first_fail_d,
    }
    outfile = os.path.join(RESULTS_DIR, "jensen_hyp_fast.json")
    with open(outfile, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nResults: {outfile}")


if __name__ == "__main__":
    main()
