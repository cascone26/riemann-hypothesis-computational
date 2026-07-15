#!/usr/bin/env python3
"""
null_space_anatomy.py — Track the near-null SUBSPACE of tau_N (not just min eigenvector).

The min eigenvector flips between two types across N due to near-degeneracy.
Fix: track the full SUBSPACE spanned by all eigenvectors with |lambda| < THRESH.

Outputs:
  1. Near-null dimension vs N
  2. Subspace-projector overlap Tr(Pi_N Pi_{N'}) / dim(Pi_N) between consecutive N
     (1.0 = subspaces identical, <1 = rotating)
  3. Symmetric vs antisymmetric decomposition of the subspace
  4. Mode content: which n-modes span the near-null space?
  5. For each near-null eigenvector: A, B, P, lambda values
"""

import json, os
import numpy as np
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
LAMBDA_SQ = 14
L  = float(mp.log(LAMBDA_SQ))
sh = float(np.sinh(L/4)**2)
THRESH = 1e-3   # near-null threshold

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/null_space_anatomy.json")


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


def build_tau(N):
    VM  = vm_table(LAMBDA_SQ)
    idx = list(range(-N, N+1))
    D   = len(idx)
    R   = TWO_PI * N / L + 80.0
    G   = 60000
    grid = np.linspace(-R, R, G)
    h    = grid[1] - grid[0]
    w    = np.full(G, h); w[0] = w[-1] = 0.5*h
    Om   = np.real(_digamma(0.25 + 0.5j * grid)) - LOGPI

    E = np.empty((D, G), dtype=complex)
    for k, n in enumerate(idx):
        s  = TWO_PI * n / L
        d  = s - grid
        sm = np.abs(d) < 1e-12
        E[k, sm]  = np.sqrt(L)
        dd         = d[~sm]
        E[k, ~sm] = (1/np.sqrt(L)) * np.exp(0.5j*dd*L) * 2*np.sin(0.5*dd*L) / dd

    arch = ((E * (w * Om)) @ E.conj().T).real / (2*np.pi)
    ns   = np.array(idx, dtype=float)
    B    = (32*L*sh * (L**2 - 16*np.pi**2*np.outer(ns,ns)) /
            np.outer(L**2+16*np.pi**2*ns**2, L**2+16*np.pi**2*ns**2))
    log_k = np.array([np.log(k) for k, _ in VM])
    wts   = np.array([lk * k**(-0.5) for k, lk in VM])
    P = np.zeros((D, D))
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]
            q = np.where(n == m,
                2*(1-np.abs(log_k)/L)*np.cos(TWO_PI*n*log_k/L),
                (np.sin(TWO_PI*m*log_k/L)-np.sin(TWO_PI*n*log_k/L))/(np.pi*(n-m)+1e-300))
            P[ii,jj] = P[jj,ii] = float(wts @ q)

    tau = arch + B - P
    return (tau+tau.T)/2, arch, B, P, idx


def subspace_projector(evec, mask):
    """Projector onto subspace spanned by eigenvectors where mask is True."""
    V = evec[:, mask]
    return V @ V.T


def subspace_overlap(Pi1, Pi2, dim1):
    """Tr(Pi1 @ Pi2) / dim1 — how much of subspace 1 lies in subspace 2."""
    return float(np.trace(Pi1 @ Pi2)) / dim1


def mode_content(evec, mask, idx, top=10):
    """For each near-null eigenvector, find dominant mode indices."""
    results = []
    cols = np.where(mask)[0]
    for c in cols:
        v = evec[:, c]
        pw = v**2
        top_idx = np.argsort(pw)[-top:][::-1]
        modes = [{'n': idx[i], 'weight': float(pw[i])} for i in top_idx if pw[i] > 0.01]
        results.append(modes)
    return results


def sym_antisym_split(evec, mask, N):
    """Project each near-null eigenvector onto sym (c[n]=c[-n]) and antisym subspaces."""
    center = N
    cols   = np.where(mask)[0]
    sym_fracs = []
    for c in cols:
        v  = evec[:, c]
        D  = len(v)
        vs = np.zeros(D); va = np.zeros(D)
        for i in range(D):
            j = D - 1 - i   # index of -idx[i]
            vs[i] = (v[i] + v[j]) / 2
            va[i] = (v[i] - v[j]) / 2
        sym_frac = float(np.dot(vs, vs))
        sym_fracs.append(sym_frac)
    return sym_fracs


def main():
    Ns   = [32, 48, 64, 96, 128, 160]
    data = []
    prev_Pi  = None
    prev_dim = None
    prev_N   = None

    print(f"Near-null subspace anatomy (threshold = {THRESH})")
    print(f"{'N':>4}  {'dim':>4}  {'lam_min':>12}  {'lam_max_nearNull':>16}  "
          f"{'overlap_prev':>13}  {'sym_fracs'}")

    for N in Ns:
        tau, arch, B, P, idx = build_tau(N)
        ev, evec = np.linalg.eigh(tau)

        mask    = np.abs(ev) < THRESH
        dim     = int(mask.sum())
        lam_min = float(ev[0])
        lam_nn  = float(ev[mask].max()) if dim > 0 else float('nan')

        Pi = subspace_projector(evec, mask)

        overlap = None
        if prev_Pi is not None and dim > 0:
            # Restrict to common (smaller) size by zero-padding the old projector
            D_new = Pi.shape[0]
            D_old = prev_Pi.shape[0]
            if D_old < D_new:
                pad = D_new - D_old
                Pi_old = np.pad(prev_Pi, ((0,pad),(0,pad)))
            else:
                Pi_old = prev_Pi[:D_new, :D_new]
            overlap = subspace_overlap(Pi_old, Pi, min(prev_dim, dim))

        mc    = mode_content(evec, mask, idx, top=5)
        sfracs = sym_antisym_split(evec, mask, N)

        # For each near-null evec: A, B, P values
        near_null_evs = []
        cols = np.where(mask)[0]
        for ci, c in enumerate(cols):
            v = evec[:, c]
            near_null_evs.append({
                'lambda': float(ev[c]),
                'A': float(v @ arch @ v),
                'B': float(v @ B @ v),
                'P': float(v @ P @ v),
                'sym_frac': sfracs[ci],
                'top_modes': mc[ci][:3]
            })

        overlap_str = f"{overlap:.4f}" if overlap is not None else "  --"
        sym_str = ' '.join(f"{s:.2f}" for s in sfracs[:5])
        print(f"{N:>4}  {dim:>4}  {lam_min:>12.4e}  {lam_nn:>16.4e}  "
              f"{overlap_str:>13}  [{sym_str}]")

        data.append({
            'N': N, 'dim': dim, 'lam_min': lam_min, 'lam_max_near_null': lam_nn,
            'subspace_overlap_with_prev': overlap,
            'near_null_evecs': near_null_evs
        })

        prev_Pi  = Pi
        prev_dim = dim
        prev_N   = N

    # Final: show the full spectrum of mode weights across near-null subspace at N=128
    print("\n=== N=128: mode content of the near-null subspace ===")
    tau, arch, B, P, idx = build_tau(128)
    ev, evec = np.linalg.eigh(tau)
    mask = np.abs(ev) < THRESH
    # aggregate mode weights across all near-null eigenvectors
    D = len(idx)
    total_wt = np.zeros(D)
    for c in np.where(mask)[0]:
        total_wt += evec[:, c]**2
    # top 20 modes by total weight
    top20 = np.argsort(total_wt)[-20:][::-1]
    print(f"  Total near-null dim = {mask.sum()}")
    for i in top20:
        n = idx[i]
        sn = TWO_PI * n / L
        print(f"    n={n:>4}  s={sn:>8.3f}  total_weight={total_wt[i]:.4f}")

    with open(OUTFILE, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"\nResults: {OUTFILE}")


if __name__ == "__main__":
    main()
