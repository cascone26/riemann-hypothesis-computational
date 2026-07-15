#!/usr/bin/env python3
"""
cvs_benchmark.py — Reproduce Groskin 2605.20224 Table 3 using connes-cvs 0.2.2.

Target: c=13, N=100, T=400, dps=80 → lam_min ≈ 2.865e-59, |gamma_1 error| ≈ 2.0e-55.

Also runs c=17, c=23 to check 113-OOM trend.
"""
import sys, json, os, time
import mpmath as mp

try:
    from connes_cvs import build_galerkin_matrix, compute_ground_state, extract_zeros
except ImportError:
    sys.exit("ERROR: connes-cvs not installed. Run: pip install connes-cvs")

GAMMA1_TRUE = mp.mpf("14.134725141734693790583409148232718309218880571040027508658")
GAMMA2_TRUE = mp.mpf("21.022039638771554992628479593896902777334340524902781754629")
GAMMA3_TRUE = mp.mpf("25.010857580145688763213790992562821818659549672557996672496")

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/cvs_benchmark.json")


def run_c(c, N, T, dps, n_zeros=5):
    mp.mp.dps = dps + 20  # extra guard digits
    L = float(mp.log(c))

    print(f"\nc={c}, N={N}, T={T}, dps={dps}, L={L:.4f}")
    t0 = time.time()
    Q = build_galerkin_matrix(c=c, N=N, T=T, dps=dps)
    t1 = time.time()
    print(f"  Q build: {t1-t0:.1f}s")

    lam_min, eigvec = compute_ground_state(Q)
    t2 = time.time()
    print(f"  eig: {t2-t1:.1f}s  lam_min={mp.nstr(lam_min, 8)}")

    # Extract zeros
    print(f"  Extracting {n_zeros} zeros...")
    zeros = extract_zeros(eigvec, L=L, n_zeros=n_zeros, dps=dps)
    t3 = time.time()
    print(f"  extract: {t3-t2:.1f}s")

    results = []
    for z in zeros:
        k        = z['k']
        gd       = z.get('gamma_detected')
        gt       = z.get('gamma_true')
        err      = z.get('error')
        err_str  = mp.nstr(err, 6) if err is not None else 'None'
        gd_str   = mp.nstr(gd, 12) if gd is not None else 'None'
        print(f"    zero {k}: detected={gd_str}  error={err_str}")
        results.append({
            'k': k,
            'gamma_true':     str(gt) if gt is not None else None,
            'gamma_detected': str(gd) if gd is not None else None,
            'error':          str(err) if err is not None else None,
        })

    sys.stdout.flush()
    return {
        'c': c, 'N': N, 'T': T, 'dps': dps, 'L': L,
        'lam_min': str(lam_min),
        'lam_min_log10': float(mp.log10(abs(lam_min))) if lam_min != 0 else None,
        'zeros': results,
        'dt_total': round(t3 - t0, 1),
    }


def main():
    all_results = []

    # Groskin Table 3 benchmark points
    cases = [
        (13, 100, 400, 80, 3),
        (17, 100, 400, 80, 3),
        (23, 100, 400, 150, 5),
    ]

    for c, N, T, dps, nz in cases:
        try:
            r = run_c(c, N, T, dps, nz)
            all_results.append(r)
        except Exception as e:
            print(f"ERROR c={c}: {e}")
            import traceback; traceback.print_exc()

    with open(OUTFILE, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()
