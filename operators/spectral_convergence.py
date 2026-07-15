#!/usr/bin/env python3
"""
spectral_convergence.py — Quantify how P and A+B spectra converge as lambda_sq grows.

Measures: ||A+B - P||_F / ||A+B||_F (relative Frobenius norm of tau)
and subspace alignment between the top-K eigenvectors of P and A+B.

This quantifies the Weil explicit formula: P → A+B as lambda_sq → infinity.
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
                       "../results/spectral_convergence.json")


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


def build_AB_P(lambda_sq):
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
    lsq_list = list(range(4, 21)) + [30, 50, 100, 200]
    results = []
    t0 = time.time()

    print(f"{'lsq':>5}  {'rel_cancel':>12}  {'P_max/AB_max':>13}  {'P_fro/AB_fro':>13}  {'dim':>5}  {'lam_min':>12}")
    for lsq in lsq_list:
        AB, P = build_AB_P(float(lsq))
        tau = AB - P

        fro_AB = np.linalg.norm(AB, 'fro')
        fro_P  = np.linalg.norm(P,  'fro')
        fro_tau= np.linalg.norm(tau,'fro')

        ev_tau = scipy.linalg.eigvalsh(tau)
        ev_P   = scipy.linalg.eigvalsh(P)
        ev_AB  = scipy.linalg.eigvalsh(AB)

        lam_min  = float(ev_tau[0])
        dim      = int((np.abs(ev_tau) < 0.02).sum())
        P_max    = float(ev_P[-1])
        AB_max   = float(ev_AB[-1])
        rel_cancel = fro_tau / fro_AB

        print(f"{lsq:>5}  {rel_cancel:>12.6f}  {P_max/AB_max:>13.6f}  {fro_P/fro_AB:>13.6f}  {dim:>5}  {lam_min:>12.4e}")
        results.append({
            'lsq': lsq, 'L': round(float(mp.log(lsq)), 6),
            'fro_tau': fro_tau, 'fro_AB': fro_AB, 'fro_P': fro_P,
            'rel_cancellation': rel_cancel,
            'P_fro_over_AB_fro': fro_P / fro_AB,
            'P_max': P_max, 'AB_max': AB_max, 'P_max_over_AB_max': P_max / AB_max,
            'dim': dim, 'lam_min': lam_min,
        })

    dt = time.time() - t0
    with open(OUTFILE, 'w') as f:
        json.dump({'N': N, 'G': G, 'results': results}, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {dt:.1f}s")


if __name__ == "__main__":
    main()
