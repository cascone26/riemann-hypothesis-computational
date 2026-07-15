#!/usr/bin/env python3
"""
dim_noprime.py — Test dim formula with P=0 (no prime term).

If dim = 2*(lambda_sq - 2) persists with P=0, the integer counting structure
is driven entirely by A (arch) + B (boundary), NOT by prime powers.
This would be a major insight into the geometric origin of the formula.

Compares:
  - tau with P included (standard)
  - tau_noP = A + B (P=0)

at lambda_sq = 4..22, step = 0.5, N=64.
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
N      = 64

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/dim_noprime.json")

THRESH = 0.02


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


def build_AB_and_P(lambda_sq, N):
    """Returns (arch+B, P) separately."""
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
    AB = (arch + B + (arch + B).T) / 2.0

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
    P = (P + P.T) / 2.0
    return AB, P


def main():
    lsq_arr = list(np.arange(4.0, 22.1, 0.5))
    results = []

    print(f"{'lsq':>6}  {'dim_full':>10}  {'dim_AB':>8}  {'lam_full':>12}  {'lam_AB':>12}")
    t0 = time.time()
    for lsq in lsq_arr:
        AB, P = build_AB_and_P(lsq, N)
        tau_full = AB - P
        tau_noP  = AB

        ev_full = scipy.linalg.eigvalsh(tau_full)
        ev_AB   = scipy.linalg.eigvalsh(tau_noP)

        dim_full = int((np.abs(ev_full) < THRESH).sum())
        dim_AB   = int((np.abs(ev_AB) < THRESH).sum())
        lam_full = float(ev_full[0])
        lam_AB   = float(ev_AB[0])
        expected = 2 * (int(lsq) - 2)

        print(f"{lsq:>6.1f}  {dim_full:>10}  {dim_AB:>8}  {lam_full:>12.4e}  {lam_AB:>12.4e}  "
              f"expect={expected}  full_ok={dim_full==expected}  AB_ok={dim_AB==expected}", flush=True)
        results.append({
            'lsq': round(lsq, 2),
            'dim_full': dim_full,
            'dim_AB': dim_AB,
            'lam_full': lam_full,
            'lam_AB': lam_AB,
            'expected': expected,
        })

    with open(OUTFILE, 'w') as f:
        json.dump({'N': N, 'G': G, 'thresh': THRESH, 'results': results}, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
