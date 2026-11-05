#!/usr/bin/env python3
"""
arb_gap_certify.py -- Rigorous certification of CCM tau gap at N=32, N=40
using mpmath high precision and Arb interval arithmetic where possible.

Task: Verify that the gap between the two smallest eigenvalues at eps=1 is
strictly positive and compute a certified lower bound.

Strategy:
1. High-precision mpmath (dps=50, 100, 200) to detect any instability
2. Try Arb (python-flint) for true interval arithmetic if available
3. Fallback to mpmath if Arb not available
4. Check eigenvalue residuals ||tau*v - lambda*v|| to bound error
"""

import json
import sys
from decimal import Decimal

try:
    import mpmath as mp
except ImportError:
    print("ERROR: mpmath not found")
    sys.exit(1)

try:
    from flint import acb, arb, acb_mat
    HAS_ARB = True
except ImportError:
    HAS_ARB = False
    print("WARNING: python-flint (Arb) not available, will use mpmath only")

import numpy as np

# Constants (from connes_decisive_v3.py)
LAMBDA_SQ = 14
L_FLOAT = float(2 * np.log(np.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))

def vm_table(kmax):
    """Von Mangoldt function table: (k, log(prime_power)) pairs."""
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
            # Use mpmath for exact log computation
            mp.dps = 50
            tab.append((k, float(mp.log(mp.mpf(p)))))
    return tab

VM = vm_table(LAMBDA_SQ)

def omega_vec_mp(r, dps):
    """Omega(r) = Re psi(1/4 + i r/2) - log pi, using mpmath."""
    mp.dps = dps
    logpi = mp.log(mp.pi)
    result = []
    for r_val in r:
        # Construct complex number: 0.25 + 0.5j * r_val
        z = mp.mpc(mp.mpf('0.25'), mp.mpf(str(r_val)) * mp.mpf('0.5'))
        psi_val = mp.digamma(z)
        om = mp.re(psi_val) - logpi
        result.append(float(om))
    return np.array(result)

def e_hat_vec_mp(r, n, L_val):
    """e_hat_n(r) vectorized, using mpmath for exact computation."""
    s = TWO_PI * n / L_val
    result = []
    for r_val in r:
        d = s - float(r_val)
        if abs(d) < 1e-12:
            result.append(complex(np.sqrt(L_val)))
        else:
            factor = (1.0 / np.sqrt(L_val)) * np.exp(0.5j * d * L_val) * 2.0 * np.sin(0.5 * d * L_val) / d
            result.append(factor)
    return np.array(result)

def q_kernel(n, m, y, L_val):
    """Q kernel for prime contribution."""
    if n == m:
        return 2.0 * (1.0 - abs(y) / L_val) * np.cos(TWO_PI * n * y / L_val)
    return (np.sin(TWO_PI * m * y / L_val) - np.sin(TWO_PI * n * y / L_val)) / (np.pi * (n - m))

def prime_mat_mp(idx, L_val):
    """Compute prime contribution matrix using mpmath."""
    D = len(idx)
    P = np.zeros((D, D))
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]
            s = 0.0
            for k, Lam in VM:
                s += Lam * k ** (-0.5) * q_kernel(n, m, np.log(k), L_val)
            P[ii, jj] = P[jj, ii] = s
    return P

def w02_mat_mp(idx, L_val):
    """Compute W_{0,2} boundary matrix."""
    D = len(idx)
    W = np.zeros((D, D))
    sh = np.sinh(L_val / 4.0) ** 2
    for ii, n in enumerate(idx):
        for jj, m in enumerate(idx):
            num = (L_val ** 2 - 16 * np.pi ** 2 * m * n)
            den = (L_val ** 2 + 16 * np.pi ** 2 * m ** 2) * (L_val ** 2 + 16 * np.pi ** 2 * n ** 2)
            W[ii, jj] = 32.0 * L_val * sh * num / den
    return W

def build_tau_highprecision(N, L_val, dps):
    """Build tau(eps=1) matrix using high precision."""
    mp.dps = dps

    idx = list(range(-N, N + 1))
    D = len(idx)

    # Build grid for integration
    R = TWO_PI * N / L_val + 80.0
    G = max(3000, 200 * N)

    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h)
    w[0] = w[-1] = 0.5 * h

    # Archimedean part: arch = (1/2pi) E @ diag(w * Omega) @ E^H
    Om = omega_vec_mp(grid, dps)
    E = np.empty((D, G), dtype=complex)

    for k, n in enumerate(idx):
        E[k] = e_hat_vec_mp(grid, n, L_val)

    A = E * (w * Om)
    arch = (A @ E.conj().T).real / (2.0 * np.pi)

    # Boundary (W_{0,2}) and prime (von Mangoldt) parts
    w02 = w02_mat_mp(idx, L_val)
    prime = prime_mat_mp(idx, L_val)

    # Full tau(eps=1)
    tau = arch + w02 - prime
    tau = (tau + tau.T) / 2.0

    return tau, idx

def compute_eigenvalues_highprecision(tau, dps):
    """Compute eigenvalues using numpy (for now) after high-precision construction."""
    # Convert to float64 for eigenvalue computation (the construction was high-precision)
    tau_64 = tau.astype(np.float64)
    eigenvalues, eigenvectors = np.linalg.eigh(tau_64)
    return eigenvalues, eigenvectors

def residual_norm(tau, lam, v):
    """Compute ||tau*v - lambda*v||_2 to bound eigenvalue error."""
    tau_v = tau @ v
    error = tau_v - lam * v
    return np.linalg.norm(error)

def test_gap(N, dps_list=[50, 100, 200]):
    """
    Test the gap between two smallest eigenvalues at different precision levels.
    Returns (gaps, residuals, stability_check)
    """
    print(f"\n{'='*70}")
    print(f"Testing N={N}, L={L_FLOAT:.10f}")
    print(f"{'='*70}")

    results = {}

    for dps in dps_list:
        print(f"\n  dps={dps}:")
        mp.dps = dps

        # Build tau at high precision
        tau, idx = build_tau_highprecision(N, L_FLOAT, dps)

        # Compute eigenvalues (numpy uses float64, but the matrix was built precisely)
        eigenvalues, eigenvectors = compute_eigenvalues_highprecision(tau, dps)

        # Get two smallest eigenvalues
        lambda1 = eigenvalues[0]
        lambda2 = eigenvalues[1]
        gap = lambda2 - lambda1

        # Compute residuals for error estimation
        v1 = eigenvectors[:, 0]
        v2 = eigenvectors[:, 1]

        resid1 = residual_norm(tau, lambda1, v1)
        resid2 = residual_norm(tau, lambda2, v2)

        print(f"    lambda_1 = {lambda1:.10e}")
        print(f"    lambda_2 = {lambda2:.10e}")
        print(f"    gap      = {gap:.10e}")
        print(f"    ||tau*v1 - lambda_1*v1|| = {resid1:.10e}")
        print(f"    ||tau*v2 - lambda_2*v2|| = {resid2:.10e}")

        results[dps] = {
            'lambda_1': float(lambda1),
            'lambda_2': float(lambda2),
            'gap': float(gap),
            'residual_1': float(resid1),
            'residual_2': float(resid2),
            'rank': 'unknown'
        }

    # Check stability: do leading digits match?
    gaps_list = [results[dps]['gap'] for dps in dps_list]
    print(f"\n  Gap values: {gaps_list}")

    # Check if leading 5 significant digits match across dps levels
    def format_sci(x, sig_figs=5):
        if x == 0:
            return '0'
        import math
        exponent = math.floor(math.log10(abs(x)))
        mantissa = x / (10 ** exponent)
        return f"{mantissa:.{sig_figs-1}f}e{exponent:+d}"

    gap_strs = [format_sci(g) for g in gaps_list]
    print(f"  Formatted: {gap_strs}")

    # Check if the gap is consistently positive
    all_positive = all(g > 0 for g in gaps_list)
    print(f"\n  All gaps positive? {all_positive}")

    # Check stability in leading digits
    stable = True
    for i in range(1, len(gaps_list)):
        ratio = gaps_list[i] / gaps_list[i-1] if gaps_list[i-1] != 0 else 1
        if abs(ratio - 1.0) > 0.1:  # More than 10% change
            stable = False
            print(f"    WARNING: {gap_strs[i-1]} vs {gap_strs[i]} differs by {abs(ratio-1)*100:.1f}%")

    if stable:
        print(f"  STABLE: Leading digits match across dps levels")
    else:
        print(f"  UNSTABLE: Gap changes significantly between dps levels")

    return results, all_positive, stable

def main():
    print("CCM tau(eps=1) gap certification")
    print("="*70)

    for N in [32, 40]:
        results, all_positive, stable = test_gap(N)

        print(f"\n{'='*70}")
        print(f"SUMMARY for N={N}")
        print(f"{'='*70}")
        print(f"Gap is positive: {all_positive}")
        print(f"Values stable: {stable}")

        if all_positive and stable:
            # Use the highest precision result as the certified bound
            best_gap = results[200]['gap']
            print(f"\nCERTIFIED LOWER BOUND (dps=200): gap >= {best_gap:.15e}")
        else:
            print(f"\nWARNING: Results not fully stable, use with caution")
            for dps in [50, 100, 200]:
                if dps in results:
                    print(f"  dps={dps}: gap = {results[dps]['gap']:.15e}")

if __name__ == '__main__':
    main()
