#!/usr/bin/env python3
"""
cvs_sweep_olympus.py — Full c=13..67 sweep reproducing Groskin Table 3.
Run on Olympus: python3 cvs_sweep_olympus.py

Output: ~/riemann/results/cvs_sweep_olympus.json
"""
import sys, json, os, time
import mpmath as mp

sys.path.insert(0, os.path.expanduser('~/riemann'))
from connes_cvs import build_galerkin_matrix, compute_ground_state, extract_zeros

GAMMA1_TRUE = mp.mpf("14.134725141734693790583409148232718309218880571040027508658")

CUTOFFS = [13, 14, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67]

OUTFILE = os.path.expanduser("~/riemann/results/cvs_sweep_olympus.json")
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)


def run_c(c, N=100, T=800, dps=150, n_zeros=5):
    mp.mp.dps = dps + 20
    L = float(mp.log(c))

    print(f"\nc={c}, N={N}, T={T}, dps={dps}, L={L:.4f}", flush=True)
    t0 = time.time()
    Q = build_galerkin_matrix(c=c, N=N, T=T, dps=dps)
    t1 = time.time()
    print(f"  Q build: {t1-t0:.1f}s", flush=True)

    lam_min, eigvec = compute_ground_state(Q)
    t2 = time.time()
    log10_lam = float(mp.log10(abs(lam_min))) if lam_min != 0 else None
    print(f"  eig: {t2-t1:.1f}s  lam_min={mp.nstr(lam_min, 8)}  log10={log10_lam:.2f}", flush=True)

    zeros = extract_zeros(eigvec, L=L, n_zeros=n_zeros, dps=dps)
    t3 = time.time()

    gamma1_error = None
    for z in zeros:
        if z['k'] == 1 and z.get('error') is not None:
            gamma1_error = str(z['error'])
            err_log10 = float(mp.log10(abs(z['error'])))
            print(f"  |gamma_1 error| = {mp.nstr(z['error'], 6)}  log10={err_log10:.2f}", flush=True)
            break

    return {
        'c': c, 'N': N, 'T': T, 'dps': dps, 'L': L,
        'lam_min': str(lam_min),
        'log10_lam_min': log10_lam,
        'gamma1_error': gamma1_error,
        'zeros': [{
            'k': z['k'],
            'gamma_true':     str(z.get('gamma_true')),
            'gamma_detected': str(z.get('gamma_detected')),
            'error':          str(z.get('error')),
        } for z in zeros],
        'dt_total': round(t3 - t0, 1),
    }


def main():
    all_results = []
    print(f"CvS sweep: {len(CUTOFFS)} cutoffs, N=100, T=800, dps=150")

    for c in CUTOFFS:
        try:
            r = run_c(c)
            all_results.append(r)
            # Save after each to preserve progress
            with open(OUTFILE, 'w') as f:
                json.dump(all_results, f, indent=2)
            print(f"  -> saved {len(all_results)}/{len(CUTOFFS)}", flush=True)
        except Exception as e:
            print(f"ERROR c={c}: {e}", flush=True)
            import traceback; traceback.print_exc()

    print(f"\nFinal save: {OUTFILE}")


if __name__ == "__main__":
    main()
