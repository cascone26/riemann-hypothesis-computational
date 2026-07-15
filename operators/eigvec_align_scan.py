#!/usr/bin/env python3
"""
eigvec_align_scan.py — Scan top-eigenvector alignment of P and A+B near lsq=12..16.

subspace_cancel.py showed top1_align jumps from ~0 at lsq=12 to 0.991 at lsq=14,
which coincides with the first Riemann zero gamma1=14.135.

This script scans lsq=11..17 at fine resolution (step=0.1) to find:
1. Exactly where the alignment transition occurs
2. Whether it coincides with gamma1=14.135 or another special value
3. Behavior of sub_cancel and other metrics through the transition
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
                       "../results/eigvec_align_scan.json")
GAMMA1 = 14.134725


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
    # Fine scan around the transition + also include gamma1 exactly
    lsq_arr = sorted(set(
        list(np.arange(11.0, 17.1, 0.1))
        + [GAMMA1, GAMMA1 - 0.5, GAMMA1 + 0.5, 12.5, 13.0, 13.5]
    ))
    results = []
    THRESH = 0.02

    print(f"{'lsq':>7}  {'dim':>4}  {'sub_cancel':>10}  {'top1':>8}  {'Pmax':>7}  {'ABmax':>7}  {'ratio':>7}")
    t0 = time.time()

    for lsq in lsq_arr:
        AB, P = build_AB_P(float(lsq))
        tau = AB - P

        ev_tau, V_tau = scipy.linalg.eigh(tau)
        dim = int((np.abs(ev_tau) < THRESH).sum())

        ev_P, V_P   = scipy.linalg.eigh(P)
        ev_AB, V_AB = scipy.linalg.eigh(AB)

        P_max  = float(ev_P[-1])
        AB_max = float(ev_AB[-1])

        # Top eigenvector alignment
        top1 = float(np.abs(V_P[:, -1] @ V_AB[:, -1]))

        # Sub-cancel in near-null
        if dim > 0:
            null_mask = np.abs(ev_tau) < THRESH
            V_null = V_tau[:, null_mask]
            AB_proj = V_null.T @ AB @ V_null
            P_proj  = V_null.T @ P  @ V_null
            fro_AB_proj = np.linalg.norm(AB_proj, 'fro')
            sub_cancel = np.linalg.norm(AB_proj - P_proj, 'fro') / fro_AB_proj
        else:
            sub_cancel = float('nan')

        marker = ' <-- GAMMA1' if abs(lsq - GAMMA1) < 0.05 else ''
        print(f"{lsq:>7.3f}  {dim:>4}  {sub_cancel:>10.6f}  {top1:>8.5f}  "
              f"{P_max:>7.4f}  {AB_max:>7.4f}  {P_max/AB_max:>7.4f}{marker}", flush=True)
        results.append({
            'lsq': round(lsq, 4),
            'dim': dim,
            'sub_cancel': float(sub_cancel),
            'top1_align': top1,
            'P_max': P_max,
            'AB_max': AB_max,
            'ratio': P_max / AB_max,
        })

    with open(OUTFILE, 'w') as f:
        json.dump({'N': N, 'G': G, 'thresh': THRESH, 'gamma1': GAMMA1, 'results': results}, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
