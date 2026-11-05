#!/usr/bin/env python3
"""
Simple test: just compute eigenvalues of a small tau(eps=1) matrix
at high precision using mpmath.
"""

import numpy as np
import mpmath as mp

# Simpler approach: just use numpy to build tau quickly, then
# use mpmath's eigsy on a converted matrix

LAMBDA_SQ = 14
L_FLOAT = float(2 * np.log(np.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))

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
            tab.append((k, float(mp.log(mp.mpf(p)))))
    return tab

VM = vm_table(LAMBDA_SQ)

def omega_vec(r):
    """Omega using scipy (fast)."""
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
    if n == m:
        return 2.0 * (1.0 - abs(y) / L_FLOAT) * np.cos(TWO_PI * n * y / L_FLOAT)
    return (np.sin(TWO_PI * m * y / L_FLOAT) - np.sin(TWO_PI * n * y / L_FLOAT)) / (np.pi * (n - m))

def prime_mat(idx):
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
    """Build tau(eps=1) using float64 (fast)."""
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

def test_eigenvalues():
    """Test eigenvalue computation at different precision levels."""
    print("Testing eigenvalue computation...")

    for N in [32, 40]:
        print(f"\n{'='*70}")
        print(f"N={N}")
        print(f"{'='*70}")

        # Build tau once using float64 (fast)
        print(f"  Building tau(eps=1) matrix...")
        tau_64, idx = build_tau(N)
        D = len(idx)
        print(f"  Matrix size: {D} x {D}")

        # Compute eigenvalues using float64 numpy first
        print(f"  Computing eigenvalues with numpy (float64)...")
        eigs_64 = np.linalg.eigvalsh(tau_64)
        print(f"    lambda_1 = {eigs_64[0]:.15e}")
        print(f"    lambda_2 = {eigs_64[1]:.15e}")
        print(f"    gap      = {eigs_64[1] - eigs_64[0]:.15e}")

        # Now try with mpmath eigenvalue solver
        for dps in [50, 100, 200]:
            print(f"\n  Converting to mpmath matrix (dps={dps})...")
            mp.dps = dps

            # Convert float64 matrix to mpmath
            tau_mp = mp.matrix([[mp.mpf(str(tau_64[i, j])) for j in range(D)] for i in range(D)])

            print(f"  Computing eigenvalues with mpmath.eigsy...")
            try:
                eigs_mp = mp.eigsy(tau_mp, eigvals_only=True)
                print(f"    lambda_1 = {mp.nstr(eigs_mp[0], 15)}")
                print(f"    lambda_2 = {mp.nstr(eigs_mp[1], 15)}")
                gap_mp = eigs_mp[1] - eigs_mp[0]
                print(f"    gap      = {mp.nstr(gap_mp, 15)}")
            except Exception as e:
                print(f"    ERROR: {e}")

if __name__ == '__main__':
    test_eigenvalues()
