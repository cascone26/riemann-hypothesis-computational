#!/usr/bin/env python3
"""
subspace_cancel.py — Study P vs A+B cancellation within the near-null subspace.

The spectral_convergence.py showed:
- P_max/AB_max → 1 (top eigenvalue converges)
- P_fro/AB_fro → 1.18+ (Frobenius DIVERGES relative to A+B)

This paradox: P and A+B don't cancel globally, but they must cancel in the near-null subspace.

This script projects P and A+B onto the near-null subspace of tau = A+B-P,
then measures:
  - ||proj(P) - proj(A+B)||_F / ||proj(A+B)||_F   (cancellation in near-null subspace)
  - Alignment of leading P eigenvectors with leading A+B eigenvectors
  - How the near-null eigenvectors of tau relate to eigenvectors of P and A+B

Tests at lsq = 4,8,12,16,20,30,50 with N=64.
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
                       "../results/subspace_cancel.json")


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
    lsq_list = [4, 6, 8, 10, 12, 14, 16, 18, 20, 30, 50]
    results = []
    THRESH = 0.02

    print(f"{'lsq':>5}  {'dim':>5}  {'sub_cancel':>12}  {'P_proj/AB_proj':>15}  "
          f"{'top1_align':>11}  {'top5_align':>11}")
    t0 = time.time()

    for lsq in lsq_list:
        AB, P = build_AB_P(float(lsq))
        tau = AB - P

        # Eigendecomposition of tau
        ev_tau, V_tau = scipy.linalg.eigh(tau)
        dim = int((np.abs(ev_tau) < THRESH).sum())

        if dim == 0:
            print(f"{lsq:>5}  {dim:>5}  (no near-null subspace)")
            continue

        # Near-null subspace: eigenvectors with |ev| < THRESH
        null_mask = np.abs(ev_tau) < THRESH
        V_null = V_tau[:, null_mask]  # (size, dim)

        # Project AB and P onto near-null subspace
        AB_proj = V_null.T @ AB @ V_null  # (dim, dim)
        P_proj  = V_null.T @ P  @ V_null  # (dim, dim)
        diff_proj = AB_proj - P_proj

        fro_AB_proj   = np.linalg.norm(AB_proj, 'fro')
        fro_diff_proj = np.linalg.norm(diff_proj, 'fro')
        sub_cancel    = fro_diff_proj / fro_AB_proj if fro_AB_proj > 0 else float('nan')
        P_proj_fro    = np.linalg.norm(P_proj, 'fro')

        # Eigenvalue alignment: how well do top eigenvectors of P and A+B agree?
        ev_P, V_P   = scipy.linalg.eigh(P)
        ev_AB, V_AB = scipy.linalg.eigh(AB)
        # Top-1 eigenvector alignment (cosine similarity)
        v_P_top  = V_P[:, -1]
        v_AB_top = V_AB[:, -1]
        top1_align = float(np.abs(v_P_top @ v_AB_top))

        # Top-5 subspace alignment (principal angles via SVD)
        K = min(5, V_P.shape[1])
        V_P5  = V_P[:, -K:]
        V_AB5 = V_AB[:, -K:]
        _, svals, _ = np.linalg.svd(V_P5.T @ V_AB5)
        top5_align = float(np.mean(svals[:K]))

        print(f"{lsq:>5}  {dim:>5}  {sub_cancel:>12.6f}  {P_proj_fro/fro_AB_proj:>15.6f}  "
              f"{top1_align:>11.6f}  {top5_align:>11.6f}")

        results.append({
            'lsq': lsq, 'dim': dim,
            'sub_cancel_rel': sub_cancel,
            'P_proj_over_AB_proj': float(P_proj_fro / fro_AB_proj),
            'top1_eigvec_align': top1_align,
            'top5_subspace_align': top5_align,
            'ev_tau_min': float(ev_tau[0]),
            'ev_P_max': float(ev_P[-1]),
            'ev_AB_max': float(ev_AB[-1]),
        })

    with open(OUTFILE, 'w') as f:
        json.dump({'N': N, 'G': G, 'thresh': THRESH, 'results': results}, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
