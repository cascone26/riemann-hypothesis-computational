#!/usr/bin/env python3
"""
p_matrix_analysis.py — Analyze the prime matrix P in isolation.

Questions:
1. What are the eigenvalues of P itself?
2. How do they compare to eigenvalues of A+B?
3. Does P have a near-null structure that mirrors dim=2*(lambda_sq-2)?

This helps explain WHY A+B ≈ P in the near-null subspace of tau = A+B-P.
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
                       "../results/p_matrix_analysis.json")


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


def build_AB_P(lambda_sq, N):
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
    AB = (arch + B + (arch + B).T) / 2.0
    P  = (P + P.T) / 2.0
    return AB, P


def main():
    results = []

    for lsq in [4.0, 6.0, 10.0, 14.0, 20.0]:
        L = float(mp.log(lsq))
        AB, P = build_AB_P(lsq, N)
        tau = AB - P

        ev_P   = scipy.linalg.eigvalsh(P)
        ev_AB  = scipy.linalg.eigvalsh(AB)
        ev_tau = scipy.linalg.eigvalsh(tau)

        # Sort all ascending
        # Spectrum stats
        n_P_small   = int((np.abs(ev_P) < 0.02).sum())
        n_AB_small  = int((np.abs(ev_AB) < 0.02).sum())
        n_tau_small = int((np.abs(ev_tau) < 0.02).sum())

        P_min, P_max = float(ev_P[0]), float(ev_P[-1])
        AB_min, AB_max = float(ev_AB[0]), float(ev_AB[-1])
        tau_min = float(ev_tau[0])

        # Projection of P eigenvectors onto AB eigenvectors (alignment)
        ev_P_top5 = [round(float(x), 6) for x in ev_P[:5]]
        ev_AB_bot5 = [round(float(x), 6) for x in ev_AB[:5]]
        ev_P_last5 = [round(float(x), 6) for x in ev_P[-5:]]
        ev_AB_last5 = [round(float(x), 6) for x in ev_AB[-5:]]

        # How closely does the MAX eigenvalue of P match MAX of AB?
        alignment = float(P_max / AB_max) if AB_max > 0 else 0.0

        print(f"\nlsq={lsq}  L={L:.4f}")
        print(f"  tau:  lam_min={tau_min:.4e}  n<0.02={n_tau_small}")
        print(f"  A+B:  lam_min={AB_min:.4e}  lam_max={AB_max:.4f}  n<0.02={n_AB_small}")
        print(f"  P:    lam_min={P_min:.4e}  lam_max={P_max:.4f}  n<0.02={n_P_small}")
        print(f"  P_max/AB_max={alignment:.6f}")
        print(f"  P top5 ev: {ev_P_last5}")
        print(f"  AB top5 ev: {ev_AB_last5}")

        results.append({
            'lsq': lsq, 'L': round(L, 6),
            'tau_min': tau_min, 'tau_dim': n_tau_small,
            'AB_min': AB_min, 'AB_max': AB_max, 'AB_dim': n_AB_small,
            'P_min': P_min, 'P_max': P_max, 'P_dim': n_P_small,
            'P_max_over_AB_max': alignment,
        })

    with open(OUTFILE, 'w') as f:
        json.dump({'N': N, 'G': G, 'results': results}, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()
