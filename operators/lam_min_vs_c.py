#!/usr/bin/env python3
"""
lam_min_vs_c.py — Study lam_min as a function of c in the CvS framework.

From the session:
  c=13 T=400/800 dps=80: lam_min ≈ 2-3e-59 (positive)
  c=23 T=400/800 dps=80: lam_min ≈ -2e-80 (at dps noise floor)
  c=50 T=400 dps=80:     lam_min ≈ -3e-27 (negative, weaker)

Questions:
1. Does lam_min → 0 as c increases (consistent with RH)?
2. What is the scaling law: lam_min ~ f(c)?
3. What is the "optimal" c for given dps and T?

Tests c = 13, 17, 23, 29, 31, 37, 41, 43 (primes) at T=400, dps=80.
"""
import sys, json, os, time
import mpmath as mp

try:
    from connes_cvs import build_galerkin_matrix, compute_ground_state
except ImportError:
    sys.exit("ERROR: connes-cvs not installed")

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/lam_min_vs_c.json")


def run_c(c, N, T, dps):
    mp.mp.dps = dps + 20
    L = float(mp.log(c))
    print(f"  c={c} L={L:.4f}...", end=' ', flush=True)
    t0 = time.time()
    Q = build_galerkin_matrix(c=c, N=N, T=T, dps=dps)
    t1 = time.time()
    lam_min, _ = compute_ground_state(Q)
    t2 = time.time()
    lam_str = mp.nstr(lam_min, 8)
    log10_lam = float(mp.log10(abs(lam_min))) if lam_min != 0 else None
    sign = '+' if lam_min >= 0 else '-'
    print(f"lam_min={sign}{mp.nstr(abs(lam_min), 4)}  "
          f"log10={log10_lam:.1f}  ({t2-t0:.0f}s)", flush=True)
    return {
        'c': c, 'N': N, 'T': T, 'dps': dps, 'L': round(L, 6),
        'lam_min': lam_str,
        'lam_min_sign': sign,
        'log10_abs_lam_min': log10_lam,
        'dt': round(t2 - t0, 1),
    }


def main():
    N = 100
    T = 400
    dps = 80

    # Test primes only (Groskin uses primes for c to avoid divisibility issues)
    c_list = [13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

    print(f"N={N}, T={T}, dps={dps}")
    print(f"{'c':>4}  {'L':>6}  {'lam_min':>16}  {'log10':>7}  {'time':>6}")

    results = []
    t0 = time.time()
    for c in c_list:
        r = run_c(c, N, T, dps)
        results.append(r)

    with open(OUTFILE, 'w') as f:
        json.dump({'N': N, 'T': T, 'dps': dps, 'results': results}, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
