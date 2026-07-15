#!/usr/bin/env python3
"""
null_basis_reconstruct.py — Reconstruct and analyze the near-null eigenvectors of tau.

The near-null subspace has dim = 2*(floor(lsq)-2)+1 vectors.
Each vector c_n (coefficients in sinc basis) defines a function:
  f(t) = Σ_n c_n * phi_n(t)  on [0, L]

This studies:
1. The structure of individual near-null eigenvectors (as functions of t)
2. Their Fourier-Mellin transforms F(z) = sqrt(L) * int_0^L f(t) e^{-izt} dt
3. Whether the near-null subspace collectively encodes ALL Riemann zeros in range
4. The evolution of eigenvectors as lambda_sq increases (how they reorganize)

Key question: Do the 2 new near-null vectors added at each integer k=floor(lsq)
correspond to the Riemann zeros in some range? Or is their structure unrelated to zeros?
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
                       "../results/null_basis_reconstruct.json")

# Riemann zeros for reference
ZEROS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
         37.586178, 40.918719, 43.327073, 48.005151, 49.773832]


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


def build_tau_eigh(lambda_sq):
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
    tau = arch + B - P
    tau = (tau + tau.T) / 2.0
    ev, V = scipy.linalg.eigh(tau)
    return ev, V, idx, float(L)


def fm_transform(c_vec, idx, L, z_grid):
    """F(z) = sqrt(L) * sum_n c_n * e^{i*alpha_n/2} * sinc(alpha_n/2) for z in z_grid."""
    alpha = TWO_PI * idx[:, None] - L * z_grid[None, :]  # (n_modes, n_z)
    sinc_part = np.where(np.abs(alpha) < 1e-10, 1.0, np.sin(alpha / 2) / (alpha / 2))
    phase_part = np.exp(0.5j * alpha)
    F = np.sqrt(L) * np.sum(c_vec[:, None] * phase_part * sinc_part, axis=0)
    return F


def main():
    # Test at several lambda_sq values
    lsq_list = [14.0, 21.0, 30.0, 50.0]
    z_grid = np.linspace(10, 55, 5000)
    THRESH = 0.02
    results = []

    t0 = time.time()
    for lsq in lsq_list:
        print(f"\nlsq={lsq}", flush=True)
        ev, V, idx, L = build_tau_eigh(lsq)
        null_mask = np.abs(ev) < THRESH
        dim = null_mask.sum()
        V_null = V[:, null_mask]  # (n_modes, dim)

        # For each near-null eigenvector, compute |F(z)| on the grid
        # Find depth of minima at known Riemann zeros
        zero_depths = []
        for j in range(min(dim, 10)):
            c = V_null[:, j]
            Fz = fm_transform(c, idx, L, z_grid)
            absF = np.abs(Fz)
            max_F = absF.max()
            # Relative depth at each known zero
            depths = {}
            for gamma in ZEROS:
                if gamma < z_grid[0] or gamma > z_grid[-1]:
                    continue
                i_near = np.argmin(np.abs(z_grid - gamma))
                # Local minimum depth
                w = 50  # window
                lo, hi = max(0, i_near-w), min(len(z_grid)-1, i_near+w)
                local_min = absF[lo:hi].min()
                depths[round(gamma, 3)] = round(float(local_min / max_F), 6)
            zero_depths.append({'eigvec_j': j, 'ev': round(float(ev[null_mask][j]), 10), 'depths': depths})

        # Also compute |F| for the MIN_KINETIC vector (weighted sum)
        ns_sq = (idx ** 2).astype(float)
        # In near-null subspace, find vector minimizing Σ n²|c_n|²
        # = find min of c^T diag(n²) c subject to c = V_null @ a, ||a||=1
        # = eigenvector of V_null^T diag(n²) V_null with smallest eigenvalue
        D = np.diag(ns_sq)
        K = V_null.T @ D @ V_null  # (dim, dim)
        ev_k, v_k = scipy.linalg.eigh(K)
        a_mink = v_k[:, 0]  # smallest kinetic eigenvector
        c_mink = V_null @ a_mink
        Fz_mink = fm_transform(c_mink, idx, L, z_grid)
        absF_mink = np.abs(Fz_mink)
        max_mink = absF_mink.max()
        mink_depths = {}
        for gamma in ZEROS:
            if gamma < z_grid[0] or gamma > z_grid[-1]:
                continue
            i_near = np.argmin(np.abs(z_grid - gamma))
            w = 50
            lo, hi = max(0, i_near-w), min(len(z_grid)-1, i_near+w)
            local_min = absF_mink[lo:hi].min()
            mink_depths[round(gamma, 3)] = round(float(local_min / max_mink), 6)

        print(f"  dim={dim}, min_kinetic zeros depths:")
        for g, d in sorted(mink_depths.items()):
            depth_str = f"{d:.4e}"
            print(f"    gamma={g:.3f}: depth={depth_str}")

        results.append({
            'lsq': lsq, 'L': round(L, 6), 'dim': int(dim),
            'min_kinetic_depths': mink_depths,
            'per_eigvec_depths': zero_depths[:5],
        })

    with open(OUTFILE, 'w') as f:
        json.dump({'N': N, 'G': G, 'zeros': ZEROS, 'results': results}, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
