#!/usr/bin/env python3
"""
near_zero_submatrix.py — Analyze the 23x23 near-zero submatrix of tau_N.

The near-zero subspace V0 = {n=0,±1,...,±11} (23 modes).
Extract tau_N restricted to V0 and analyze:
  1. Is it tridiagonal (or close)?
  2. Are all leading principal minors positive?
  3. What are its eigenvalues?
  4. What is its conditioning?
"""

import json, os, numpy as np
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))
LAMBDA_SQ = 14
L = float(mp.log(LAMBDA_SQ))
sh = float(np.sinh(L/4)**2)

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/near_zero_submatrix.json")

NZ_HALF = 11  # near-zero subspace: n=-11,...,0,...,+11 = 23 modes


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


def build_full_tau(N):
    VM = vm_table(LAMBDA_SQ)
    idx = list(range(-N, N+1))
    D = len(idx)
    R = TWO_PI * N / L + 80.0
    G = 60000
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h); w[0] = w[-1] = 0.5*h
    Om = np.real(_digamma(0.25 + 0.5j * grid)) - LOGPI

    E = np.empty((D, G), dtype=complex)
    for k, n in enumerate(idx):
        s = TWO_PI * n / L
        d = s - grid
        sm = np.abs(d) < 1e-12
        E[k, sm] = np.sqrt(L)
        dd = d[~sm]
        E[k, ~sm] = (1/np.sqrt(L)) * np.exp(0.5j*dd*L) * 2*np.sin(0.5*dd*L) / dd

    arch = ((E * (w * Om)) @ E.conj().T).real / (2*np.pi)
    ns = np.array(idx, dtype=float)
    B = 32*L*sh*(L**2 - 16*np.pi**2*np.outer(ns,ns)) / np.outer(L**2+16*np.pi**2*ns**2, L**2+16*np.pi**2*ns**2)
    log_k = np.array([np.log(k) for k, _ in VM])
    wts = np.array([lk * k**(-0.5) for k, lk in VM])
    P = np.zeros((D, D))
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]
            q = np.where(n == m,
                2*(1-np.abs(log_k)/L)*np.cos(TWO_PI*n*log_k/L),
                (np.sin(TWO_PI*m*log_k/L)-np.sin(TWO_PI*n*log_k/L))/(np.pi*(n-m)+1e-300))
            P[ii,jj] = P[jj,ii] = float(wts @ q)

    tau = arch + B - P
    return (tau + tau.T)/2, arch, B, P, idx


def analyze_submatrix(tau, arch, B, P, idx, N, nz_half=NZ_HALF):
    # Indices of near-zero subspace in the full matrix
    nz_modes = list(range(-nz_half, nz_half+1))  # -11,...,0,...,11
    nz_indices = [idx.index(n) for n in nz_modes]
    D_nz = len(nz_indices)

    # Extract submatrices
    tau_nz = tau[np.ix_(nz_indices, nz_indices)]
    arch_nz = arch[np.ix_(nz_indices, nz_indices)]
    B_nz = B[np.ix_(nz_indices, nz_indices)]
    P_nz = P[np.ix_(nz_indices, nz_indices)]

    print(f"\n=== Near-zero submatrix analysis (nz_half={nz_half}, D_nz={D_nz}) ===")
    print(f"  Modes: n={-nz_half},...,0,...,{nz_half}")

    # Diagonal
    diag = np.diag(tau_nz)
    print(f"\n  Diagonal tau_nz[i,i] (ordered by n):")
    for i, n in enumerate(nz_modes):
        print(f"    n={n:3d}: tau={diag[i]:.6f}  (arch={np.diag(arch_nz)[i]:.4f}, B={np.diag(B_nz)[i]:.4f}, P={np.diag(P_nz)[i]:.4f})")

    # Tridiagonality check: how much mass is on the diagonal vs off-diagonal bands?
    total_mass = np.sum(tau_nz**2)
    diag_mass = np.sum(np.diag(tau_nz)**2)
    super_diag_mass = np.sum(np.diag(tau_nz, 1)**2)
    sub_diag_mass = np.sum(np.diag(tau_nz, -1)**2)
    off_rest = total_mass - diag_mass - super_diag_mass - sub_diag_mass

    print(f"\n  Frobenius mass distribution:")
    print(f"    Diagonal:       {diag_mass/total_mass*100:.1f}%  (mass={diag_mass:.4f})")
    print(f"    Super-diagonal: {super_diag_mass/total_mass*100:.1f}%  (mass={super_diag_mass:.4f})")
    print(f"    Sub-diagonal:   {sub_diag_mass/total_mass*100:.1f}%  (mass={sub_diag_mass:.4f})")
    print(f"    Rest (off-tridiag): {off_rest/total_mass*100:.1f}%  (mass={off_rest:.4f})")

    # Off-diagonal elements |tau_nz[i,j]| for various distances
    print(f"\n  Off-diagonal decay |tau_nz[i,j]| vs |i-j|:")
    for dist in range(1, min(12, D_nz//2)):
        vals = [abs(tau_nz[i, i+dist]) for i in range(D_nz-dist)]
        print(f"    |i-j|={dist}: max={max(vals):.4e}, mean={np.mean(vals):.4e}")

    # Eigenvalues
    ev = np.linalg.eigvalsh(tau_nz)
    print(f"\n  Eigenvalues of tau_nz (23x23):")
    print(f"    All eigenvalues: {ev}")
    print(f"    Min: {ev[0]:.4e}, Max: {ev[-1]:.4f}")
    print(f"    Positive? {np.all(ev > 0)}")
    print(f"    All > 0 up to threshold 1e-4? {np.all(ev > -1e-4)}")

    # Leading principal minors
    print(f"\n  Leading principal minors (Sylvester criterion):")
    all_pos = True
    minors = []
    for k in range(1, D_nz+1):
        m = np.linalg.det(tau_nz[:k, :k])
        minors.append(float(m))
        pos = m > 0
        if not pos:
            all_pos = False
        if k <= 5 or k >= D_nz-2 or not pos:
            print(f"    M_{k}: {m:.6e}  {'✓' if pos else '✗ NEGATIVE'}")

    print(f"\n  All leading principal minors > 0? {all_pos}")

    # Full matrix eigenvalues for comparison
    ev_full = np.linalg.eigvalsh(tau)
    nz_full = np.sum(np.abs(ev_full) < 1e-4)
    first_gap = float(ev_full[ev_full > 1e-4][0]) if np.any(ev_full > 1e-4) else None

    return {
        'N': N, 'nz_half': nz_half, 'D_nz': D_nz,
        'tau_nz_diagonal': diag.tolist(),
        'tau_nz_eigenvalues': ev.tolist(),
        'tau_nz_min_ev': float(ev[0]), 'tau_nz_max_ev': float(ev[-1]),
        'all_ev_positive': bool(np.all(ev > 0)),
        'all_minors_positive': all_pos,
        'diag_frac': float(diag_mass/total_mass),
        'tridiag_frac': float((diag_mass+super_diag_mass+sub_diag_mass)/total_mass),
        'leading_minors': minors,
        'full_lambda_min': float(ev_full[0]),
        'full_nz_count': int(nz_full),
        'full_first_gap': first_gap
    }


def main():
    all_results = []
    for N in [32, 48, 64]:
        print(f"\nBuilding N={N}...", flush=True)
        tau, arch, B, P, idx = build_full_tau(N)
        r = analyze_submatrix(tau, arch, B, P, idx, N, nz_half=NZ_HALF)
        all_results.append(r)

    with open(OUTFILE, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults: {OUTFILE}", flush=True)


if __name__ == "__main__":
    main()
