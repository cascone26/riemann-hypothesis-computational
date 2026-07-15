#!/usr/bin/env python3
"""
dim_scaling.py — Test dim formula at large lambda_sq with varying N.

Formula hypothesis: dim(tau, thresh=0.02) = 2*(floor(lsq)-2) for integer lsq.
But N=64 may be insufficient for large lsq. This tests N=32,48,64,96,128 at lsq=20,30,50.

Key question: does dim grow linearly with lsq even for large lsq when N is large enough?
"""
import json, os, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 50
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/dim_scaling.json")


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


def build_tau(lambda_sq, N):
    L  = float(mp.log(lambda_sq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lambda_sq)
    idx = np.arange(-N, N + 1, dtype=float)
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
        lkpp = np.array([np.log(k) for k, _ in VM])
        wts  = np.array([lk * k ** (-0.5) for k, lk in VM])
        ph   = TWO_PI * np.outer(ns, lkpp) / L
        S    = np.sin(ph) @ wts
        dn   = ns[:, None] - ns[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            P = (S[None, :] - S[:, None]) / (np.pi * dn)
        P[dn == 0] = 0.0
        np.fill_diagonal(P, np.cos(ph) @ (2 * wts * (1 - lkpp / L)))
    else:
        P = np.zeros((len(idx), len(idx)))
    tau = (arch + B - P)
    tau = (tau + tau.T) / 2.0
    return tau


def main():
    N_values  = [32, 48, 64, 96, 128]
    lsq_values = [20, 30, 50, 100]
    results = []

    # Expected dim from formula 2*(lsq-2) (integer lsq)
    print(f"{'lsq':>5}  {'N':>4}  {'dim':>5}  {'exp':>5}  {'lam_min':>12}  {'dim/N2':>8}")
    t0 = time.time()
    for lsq in lsq_values:
        exp_dim = 2 * (lsq - 2)
        for N in N_values:
            t1 = time.time()
            tau = build_tau(float(lsq), N)
            ev = scipy.linalg.eigvalsh(tau)
            dim = int((np.abs(ev) < 0.02).sum())
            lam_min = float(ev[0])
            matrix_size = 2 * N + 1
            ratio = dim / matrix_size
            print(f"{lsq:>5}  {N:>4}  {dim:>5}  {exp_dim:>5}  {lam_min:>12.4e}  {ratio:>8.4f}  "
                  f"({time.time()-t1:.1f}s)", flush=True)
            results.append({
                'lsq': lsq, 'N': N, 'matrix_size': matrix_size,
                'dim': dim, 'expected': exp_dim,
                'lam_min': lam_min, 'dim_fraction': ratio,
            })

    with open(OUTFILE, 'w') as f:
        json.dump({'results': results}, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
