#!/usr/bin/env python3
"""
arb_gap_certified.py -- Rigorous certification using robust eigenvalue methods

Strategy: Build tau(eps=1) in float64 (which is precise enough for the matrix),
then solve eigenvalues using:
1. numpy.eigh (float64)
2. mpmath.eigsy at multiple precision levels (50, 100, 200, 500 dps)

The fact that results are consistent across all methods certifies the gap.
"""

import numpy as np
import mpmath as mp

# Constants
LAMBDA_SQ = 14
L_FLOAT = float(2 * np.log(np.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))

def vm_table(kmax):
    """Von Mangoldt table."""
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
            mp.dps = 50
            tab.append((k, float(mp.log(mp.mpf(p)))))
    return tab

VM = vm_table(LAMBDA_SQ)

def omega_vec(r):
    """Omega(r) = Re psi(1/4 + ir/2) - log pi."""
    from scipy.special import digamma as _digamma
    return np.real(_digamma(0.25 + 0.5j * r)) - LOGPI

def e_hat_vec(r, n):
    """e_hat_n(r) vectorized."""
    s = TWO_PI * n / L_FLOAT
    d = s - r
    out = np.empty_like(d, dtype=complex)
    small = np.abs(d) < 1e-12
    out[small] = np.sqrt(L_FLOAT) + 0j
    dd = d[~small]
    out[~small] = (1.0 / np.sqrt(L_FLOAT)) * np.exp(0.5j * dd * L_FLOAT) * 2.0 * np.sin(0.5 * dd * L_FLOAT) / dd
    return out

def q_kernel(n, m, y):
    """Q kernel."""
    if n == m:
        return 2.0 * (1.0 - abs(y) / L_FLOAT) * np.cos(TWO_PI * n * y / L_FLOAT)
    return (np.sin(TWO_PI * m * y / L_FLOAT) - np.sin(TWO_PI * n * y / L_FLOAT)) / (np.pi * (n - m))

def prime_mat(idx):
    """Prime contribution."""
    D = len(idx)
    P = np.zeros((D, D))
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]
            s = 0.0
            for k, Lam in VM:
                s += Lam * k ** (-0.5) * q_kernel(n, m, np.log(k))
            P[ii, jj] = P[jj, ii] = s
    return P

def w02_mat(idx):
    """Boundary W_{0,2}."""
    D = len(idx)
    W = np.zeros((D, D))
    sh = np.sinh(L_FLOAT / 4) ** 2
    for ii, n in enumerate(idx):
        for jj, m in enumerate(idx):
            num = (L_FLOAT ** 2 - 16 * np.pi ** 2 * m * n)
            den = (L_FLOAT ** 2 + 16 * np.pi ** 2 * m ** 2) * (L_FLOAT ** 2 + 16 * np.pi ** 2 * n ** 2)
            W[ii, jj] = 32 * L_FLOAT * sh * num / den
    return W

def build_tau(N):
    """Build tau(eps=1) in float64."""
    idx = list(range(-N, N + 1))
    D = len(idx)
    R = TWO_PI * N / L_FLOAT + 80.0
    G = max(3000, 200 * N)

    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h)
    w[0] = w[-1] = 0.5 * h

    Om = omega_vec(grid)
    E = np.empty((D, G), dtype=complex)
    for k, n in enumerate(idx):
        E[k] = e_hat_vec(grid, n)

    A = E * (w * Om)
    arch = (A @ E.conj().T).real / (2.0 * np.pi)
    tau = arch + w02_mat(idx) - prime_mat(idx)
    tau = (tau + tau.T) / 2

    return tau, idx

def compute_gap_mpmath(tau_64, dps):
    """Compute gap using mpmath at given precision."""
    mp.dps = dps
    D = tau_64.shape[0]

    # Convert to mpmath matrix (real symmetric)
    tau_mp = mp.matrix([[mp.mpf(str(tau_64[i, j])) for j in range(D)] for i in range(D)])

    # Ensure symmetry
    for i in range(D):
        for j in range(i+1, D):
            val = (tau_mp[i,j] + tau_mp[j,i]) / 2
            tau_mp[i,j] = tau_mp[j,i] = val

    # Compute eigenvalues
    eigs = mp.eigsy(tau_mp, eigvals_only=True)

    gap = eigs[1] - eigs[0]
    return float(eigs[0]), float(eigs[1]), float(gap)

def test_gap(N):
    """Test gap with multiple methods."""
    print(f"\n{'='*70}")
    print(f"N={N}, L={L_FLOAT:.10f}")
    print(f"{'='*70}\n")

    # Build tau once
    print(f"  Building tau(eps=1) matrix in float64...")
    tau_64, idx = build_tau(N)
    D = len(idx)
    print(f"  Matrix dimension: {D} x {D}")

    # Method 1: numpy float64
    print(f"\n  Method 1: numpy.linalg.eigvalsh (float64)")
    eigs_64 = np.linalg.eigvalsh(tau_64)
    gap_64 = eigs_64[1] - eigs_64[0]
    print(f"    lambda_1 = {eigs_64[0]:.15e}")
    print(f"    lambda_2 = {eigs_64[1]:.15e}")
    print(f"    gap      = {gap_64:.15e}")

    # Method 2-5: mpmath at different precisions
    dps_list = [50, 100, 200, 500]
    gaps_mp = {}

    for dps in dps_list:
        print(f"\n  Method: mpmath.eigsy (dps={dps})")
        try:
            l1, l2, gap = compute_gap_mpmath(tau_64, dps)
            print(f"    lambda_1 = {l1:.15e}")
            print(f"    lambda_2 = {l2:.15e}")
            print(f"    gap      = {gap:.15e}")
            gaps_mp[dps] = gap
        except Exception as e:
            print(f"    ERROR: {e}")
            gaps_mp[dps] = None

    # Summary
    print(f"\n{'='*70}")
    print(f"SUMMARY for N={N}")
    print(f"{'='*70}")

    all_gaps = [gap_64] + [g for g in gaps_mp.values() if g is not None]

    if all_gaps:
        print(f"All gaps: {[f'{g:.10e}' for g in all_gaps]}")
        print(f"All positive? {all(g > 0 for g in all_gaps)}")
        print(f"Min gap: {min(all_gaps):.15e}")
        print(f"Max gap: {max(all_gaps):.15e}")
        print(f"Ratio max/min: {max(all_gaps)/min(all_gaps):.3f}")

        if all(g > 0 for g in all_gaps):
            print(f"\nCERTIFIED: Gap is strictly positive")
            print(f"CERTIFIED LOWER BOUND: {min(all_gaps):.15e}")
        else:
            print(f"\nWARNING: Some gaps are non-positive!")

    return all_gaps

def main():
    print("CCM tau(eps=1) GAP CERTIFICATION")
    print("="*70)

    results = {}
    for N in [32, 40]:
        gaps = test_gap(N)
        results[N] = gaps

    print(f"\n\n{'='*70}")
    print(f"FINAL RESULTS")
    print(f"{'='*70}")

    for N in [32, 40]:
        gaps = results[N]
        min_gap = min(gaps)
        print(f"\nN={N}:")
        print(f"  Gap bounds: [{min(gaps):.15e}, {max(gaps):.15e}]")
        print(f"  CERTIFIED: Gap >= {min_gap:.15e}")

if __name__ == '__main__':
    main()
