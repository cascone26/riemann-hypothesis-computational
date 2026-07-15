#!/usr/bin/env python3
"""
large_n_scan.py — Push to large N using vectorized builder.
Key question: does the 24-dim near-null cluster stay at 24 for N >> 128?
Tracks: dim, spectral_gap, lam_min, lam_max_near_null, build_time.

Run on Olympus: python3 ~/riemann/large_n_scan.py
"""
import sys, json, time, os
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
LAMBDA_SQ = 14
THRESH = 1e-3   # near-null threshold
OUTFILE = os.path.expanduser("~/riemann/results/large_n_scan.json")
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)


def vm_table(kmax):
    tab = []
    for k in range(2, kmax + 1):
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


def build_tau_fast(N, lambda_sq=14):
    L  = float(mp.log(lambda_sq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lambda_sq)
    idx = np.arange(-N, N + 1, dtype=float)
    D   = len(idx)
    R   = TWO_PI * N / L + 80.0
    G   = 60000   # 60K grid — sufficient for near-null dim, matches original loop scripts

    grid = np.linspace(-R, R, G)
    h    = grid[1] - grid[0]
    w    = np.full(G, h); w[0] = w[-1] = 0.5 * h
    Om   = np.real(_digamma(0.25 + 0.5j * grid)) - LOGPI

    svals = TWO_PI * idx / L
    d     = svals[:, None] - grid[None, :]
    tiny  = np.abs(d) < 1e-12
    E = np.where(tiny, np.sqrt(L),
                 (1.0 / np.sqrt(L)) * np.exp(0.5j * d * L)
                 * 2.0 * np.sin(0.5 * d * L) / (d + 1e-300 * tiny))
    arch = ((E * (w * Om)) @ E.conj().T).real / (2.0 * np.pi)

    ns = idx
    B  = (32.0 * L * sh * (L**2 - 16.0 * np.pi**2 * np.outer(ns, ns)) /
          np.outer(L**2 + 16.0 * np.pi**2 * ns**2,
                   L**2 + 16.0 * np.pi**2 * ns**2))

    log_k_pp = np.array([np.log(k) for k, _ in VM])  # log(prime_power) for phases
    wts      = np.array([lk * k**(-0.5) for k, lk in VM])
    phases = TWO_PI * np.outer(ns, log_k_pp) / L
    S = np.sin(phases) @ wts
    dn = ns[:, None] - ns[None, :]   # dn[i,j] = ns[i]-ns[j]
    with np.errstate(divide='ignore', invalid='ignore'):
        P = (S[None, :] - S[:, None]) / (np.pi * dn)
    P[dn == 0] = 0.0
    diag_wts = 2.0 * wts * (1.0 - log_k_pp / L)
    np.fill_diagonal(P, np.cos(phases) @ diag_wts)

    tau = arch + B - P
    return (tau + tau.T) / 2.0


def run(N):
    t0 = time.time()
    tau = build_tau_fast(N, LAMBDA_SQ)
    t_build = time.time() - t0

    t1 = time.time()
    ev = scipy.linalg.eigvalsh(tau)
    t_eig = time.time() - t1

    mask    = np.abs(ev) < THRESH
    dim     = int(mask.sum())
    lam_min = float(ev[0])
    lam_nn  = float(ev[mask].max()) if dim > 0 else float('nan')
    # first gap = distance from top of near-null cluster to first non-near-null
    above   = ev[~mask & (ev > 0)]
    gap     = float(above.min() - ev[mask].max()) if (dim > 0 and len(above) > 0) else float('nan')

    t_total = time.time() - t0
    return {
        'N': N, 'D': 2*N+1, 'dim': dim,
        'lam_min': lam_min, 'lam_max_near_null': lam_nn,
        'spectral_gap': gap,
        't_build': round(t_build, 2), 't_eig': round(t_eig, 2), 't_total': round(t_total, 2)
    }


def main():
    Ns = [100, 150, 200, 300, 400, 500, 700, 1000]
    print(f"Large-N near-null scan | lambda_sq={LAMBDA_SQ} | thresh={THRESH}")
    print(f"{'N':>5}  {'D':>5}  {'dim':>4}  {'lam_min':>12}  {'lam_max_nn':>12}  "
          f"{'gap':>10}  {'t_build':>8}  {'t_eig':>8}")
    results = []
    for N in Ns:
        try:
            r = run(N)
            print(f"{r['N']:>5}  {r['D']:>5}  {r['dim']:>4}  {r['lam_min']:>12.4e}  "
                  f"{r['lam_max_near_null']:>12.4e}  {r['spectral_gap']:>10.4e}  "
                  f"{r['t_build']:>8.2f}s  {r['t_eig']:>8.2f}s")
            sys.stdout.flush()
            results.append(r)
        except Exception as e:
            print(f"N={N} FAILED: {e}")
            sys.stdout.flush()

    with open(OUTFILE, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()
