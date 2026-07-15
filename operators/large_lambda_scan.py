#!/usr/bin/env python3
"""
large_lambda_scan.py — Push lambda^2 to large values on Olympus.

Tests: does tau remain PSD (lam_min ~ 0) as lambda^2 grows?
At large lambda^2, hundreds of zeros are "inside" the window.
Key question: does the spectral gap stay open? Does dim follow any pattern?

Run on Olympus: python3 ~/riemann/large_lambda_scan.py
"""
import sys, json, time, os
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 50
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000
N      = 100     # moderate N — fast, reliable gap detection
THRESH = 1e-2    # wider threshold since gap may shrink at large lambda^2
OUTFILE = os.path.expanduser("~/riemann/results/large_lambda_scan.json")
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)


def vm_table(kmax):
    tab = []
    for k in range(2, int(kmax) + 1):
        n, p, f, ok = k, 0, 2, True
        while f * f <= n:
            if n % f == 0:
                if p == 0: p = f
                elif p != f: ok = False; break
                while n % f == 0: n //= f
            f += 1
        if n > 1:
            if p == 0: p = n
            elif p != n: ok = False
        if ok and p > 0:
            tab.append((k, float(mp.log(p))))
    return tab


def build_tau(lambda_sq, N, G=60000):
    L  = float(mp.log(lambda_sq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lambda_sq)
    idx = np.arange(-N, N + 1, dtype=float)
    D   = len(idx)
    ns  = idx
    R   = TWO_PI * N / L + 80.0
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]; w = np.full(G, h); w[0] = w[-1] = 0.5 * h
    Om = np.real(_digamma(0.25 + 0.5j * grid)) - LOGPI
    d = (TWO_PI * ns / L)[:, None] - grid[None, :]
    tiny = np.abs(d) < 1e-12
    E = np.where(tiny, np.sqrt(L),
                 (1 / np.sqrt(L)) * np.exp(0.5j * d * L)
                 * 2 * np.sin(0.5 * d * L) / (d + 1e-300 * tiny))
    arch = ((E * (w * Om)) @ E.conj().T).real / (2.0 * np.pi)
    B = (32.0 * L * sh * (L ** 2 - 16.0 * np.pi ** 2 * np.outer(ns, ns)) /
         np.outer(L ** 2 + 16.0 * np.pi ** 2 * ns ** 2,
                  L ** 2 + 16.0 * np.pi ** 2 * ns ** 2))
    if VM:
        log_k_pp = np.array([np.log(k) for k, _ in VM])
        wts      = np.array([lk * k ** (-0.5) for k, lk in VM])
        phases   = TWO_PI * np.outer(ns, log_k_pp) / L
        S        = np.sin(phases) @ wts
        dn       = ns[:, None] - ns[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            P = (S[None, :] - S[:, None]) / (np.pi * dn)
        P[dn == 0] = 0.0
        np.fill_diagonal(P, np.cos(phases) @ (2 * wts * (1 - log_k_pp / L)))
    else:
        P = np.zeros((D, D))
    tau = arch + B - P
    return (tau + tau.T) / 2.0


# Approximate zero counts: N(T) ~ T/(2*pi) * log(T/(2*pi*e))
def approx_zero_count(T):
    if T <= 0: return 0
    return T / (2 * np.pi) * np.log(T / (2 * np.pi * np.e))


def main():
    # lambda_sq values: pick values that span many zeros
    # N(T) zeros up to height T:
    #   T=14.1: ~0 zeros  (just below first)
    #   T=50:   ~11 zeros
    #   T=100:  ~29 zeros
    #   T=500:  ~201 zeros
    #   T=1000: ~421 zeros
    #   T=5000: ~2396 zeros
    lambda_sqs = [
        16, 20, 25, 36, 50, 64,
        100, 150, 200, 300, 500,
        750, 1000, 2000, 5000, 10000,
    ]

    results = []
    print(f"Large-lambda scan | N={N} | G={G} | thresh={THRESH}")
    print(f"{'lam_sq':>8}  {'n_pp':>5}  {'~n_zeros':>8}  {'dim':>5}  "
          f"{'lam_min':>14}  {'gap':>12}  {'t':>6}")

    for lsq in lambda_sqs:
        t0 = time.time()
        try:
            tau = build_tau(lsq, N, G)
            ev  = scipy.linalg.eigvalsh(tau)
            dt  = time.time() - t0
            mask = np.abs(ev) < THRESH
            dim  = int(mask.sum())
            lam_min = float(ev[0])
            above = ev[~mask & (ev > 0)]
            gap = float(above.min() - ev[mask].max()) if (dim > 0 and len(above) > 0) else float('nan')
            n_pp = len(vm_table(lsq))
            n_zeros_approx = int(approx_zero_count(lsq))

            print(f"{lsq:>8}  {n_pp:>5}  {n_zeros_approx:>8}  {dim:>5}  "
                  f"{lam_min:>14.4e}  {gap:>12.4e}  {dt:>6.2f}s")
            sys.stdout.flush()
            results.append({
                'lambda_sq': lsq, 'N': N, 'n_pp': n_pp,
                'n_zeros_approx': n_zeros_approx,
                'dim': dim, 'lam_min': lam_min, 'gap': gap,
                't': round(dt, 2)
            })
        except Exception as e:
            print(f"lam_sq={lsq} ERROR: {e}")
            sys.stdout.flush()

    with open(OUTFILE, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()
