#!/usr/bin/env python3
"""
gershgorin_analysis.py — Row sum (Gershgorin) structure of tau_N.

For each row i: compute R_i = Σ_{j≠i} |tau[i,j]|  (Gershgorin radius)
Gershgorin: all eigenvalues in union of disks [tau[i,i] - R_i, tau[i,i] + R_i]

Key questions:
  1. Is tau[i,i] - R_i > 0 for any rows? (diagonal dominance → eigenvalue lower bound)
  2. How tight is the bound? How does gap tau[i,i] - R_i vary with i and N?
  3. Which rows have the smallest gap (i.e., which modes are most "at risk")?
  4. What's the dominant contribution to R_i? (arch, B, or P off-diagonal?)
"""

import json, os, time
import numpy as np
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))
LAMBDA_SQ = 14
L = float(mp.log(LAMBDA_SQ))
sh = float(np.sinh(L/4)**2)

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/gershgorin_analysis.json")


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


def build_tau_components(N):
    """Build tau_N and return arch, B, P separately."""
    VM = vm_table(LAMBDA_SQ)
    idx = list(range(-N, N+1))
    D = len(idx)
    R = TWO_PI * N / L + 80.0
    G = max(60_000, 3000 * N)
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

    tau = (arch + B - P + (arch + B - P).T) / 2
    return tau, arch, B, P


def analyze_gershgorin(tau, arch, B, P, N, top_k=20):
    """Analyze Gershgorin structure row by row."""
    D = tau.shape[0]
    idx = list(range(-N, N+1))

    # Per-row analysis
    diag = np.diag(tau)
    off_sums = np.array([np.sum(np.abs(tau[i])) - np.abs(tau[i,i]) for i in range(D)])
    gaps = diag - off_sums  # positive gap = Gershgorin row dominant

    # Same for components
    diag_arch = np.diag(arch)
    off_arch = np.array([np.sum(np.abs(arch[i])) - np.abs(arch[i,i]) for i in range(D)])
    diag_B = np.diag(B)
    off_B = np.array([np.sum(np.abs(B[i])) - np.abs(B[i,i]) for i in range(D)])
    diag_P = np.diag(P)
    off_P = np.array([np.sum(np.abs(P[i])) - np.abs(P[i,i]) for i in range(D)])

    # Eigenvalues for reference
    ev = np.linalg.eigvalsh(tau)

    print(f"\n=== Gershgorin Analysis N={N} ===")
    print(f"  lambda_min={ev[0]:.4e}, lambda_max={ev[-1]:.4f}")
    print(f"  tau[0,0] = {diag[N]:.6f}, off_sum[0] = {off_sums[N]:.6f}, gap[0] = {gaps[N]:.6f}")
    print(f"  Rows with gap > 0: {np.sum(gaps > 0)} / {D}")
    print(f"  Min gap (Gershgorin lower bound): {np.min(gaps):.4e}")
    print(f"  Max gap: {np.max(gaps):.4f}")
    print()
    print(f"  {'n':>5}  {'diag':>9}  {'off_sum':>9}  {'gap':>9}  {'arch_off':>9}  {'B_off':>9}  {'P_off':>9}")

    # Show worst rows (smallest gap)
    worst_idx = np.argsort(gaps)[:top_k]
    for i in worst_idx:
        n = idx[i]
        print(f"  {n:>5}  {diag[i]:>9.5f}  {off_sums[i]:>9.5f}  {gaps[i]:>9.4e}  "
              f"{off_arch[i]:>9.5f}  {off_B[i]:>9.5f}  {off_P[i]:>9.5f}")

    # The n=0 row specifically
    i0 = N  # center index
    print(f"\n  Row n=0 decomposition:")
    print(f"    tau[0,0] = arch[0,0] + B[0,0] - P[0,0]")
    print(f"             = {diag_arch[i0]:.6f} + {diag_B[i0]:.6f} - {diag_P[i0]:.6f} = {diag[i0]:.6f}")
    print(f"    off_sum(tau[0,j]) = {off_sums[i0]:.6f}")
    print(f"      from arch: {off_arch[i0]:.6f}")
    print(f"      from B:    {off_B[i0]:.6f}")
    print(f"      from P:    {off_P[i0]:.6f}")
    print(f"    Gap (diagonal dominance margin): {gaps[i0]:.6e}")

    # Check if diagonally dominant for any rows
    n_dominant = np.sum(gaps > 0)
    if n_dominant > 0:
        print(f"\n  {n_dominant} rows are diagonally dominant (Gershgorin bound > 0)")
    else:
        print(f"\n  NO rows are diagonally dominant — all gaps < 0")
        print(f"  Min gap = {np.min(gaps):.4e} at n={idx[np.argmin(gaps)]}")

    return {
        'N': N,
        'lambda_min': float(ev[0]),
        'n_dominant': int(n_dominant),
        'min_gap': float(np.min(gaps)),
        'max_gap': float(np.max(gaps)),
        'gap_n0': float(gaps[i0]),
        'diag_n0': float(diag[i0]),
        'off_sum_n0': float(off_sums[i0]),
        'worst_rows': [{'n': idx[i], 'diag': float(diag[i]), 'off_sum': float(off_sums[i]),
                        'gap': float(gaps[i])} for i in worst_idx[:10]],
        'gaps_all': gaps.tolist()
    }


def main():
    all_results = []

    for N in [16, 32, 48, 64]:
        print(f"\nBuilding N={N}...", flush=True)
        t0 = time.time()
        tau, arch, B, P = build_tau_components(N)
        t1 = time.time()
        print(f"  Built in {t1-t0:.1f}s", flush=True)
        r = analyze_gershgorin(tau, arch, B, P, N)
        all_results.append(r)

        with open(OUTFILE, 'w') as f:
            json.dump(all_results, f, indent=2)

    print(f"\nResults: {OUTFILE}", flush=True)


if __name__ == "__main__":
    main()
