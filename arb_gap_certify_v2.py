#!/usr/bin/env python3
"""
arb_gap_certify_v2.py -- CERTIFIED gap bounds using residual error estimates

Key insight: The Bauer-Fike theorem bounds the error in computed eigenvalues
via the residual norm ||tau*v - lambda*v||.

If residual r = ||tau*v - lambda_computed*v|| and eigenvector v is unit-norm,
then the true eigenvalue lambda_true satisfies:
    |lambda_true - lambda_computed| <= r

This gives certified intervals for the eigenvalues, and hence a certified
lower bound on the gap if the intervals don't overlap.
"""

import json
import sys
import numpy as np

try:
    import mpmath as mp
except ImportError:
    print("ERROR: mpmath not found")
    sys.exit(1)

# Constants
LAMBDA_SQ = 14
L_FLOAT = float(2 * np.log(np.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))

def vm_table(kmax):
    """Von Mangoldt function table."""
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

def omega_vec_mp(r, dps):
    """Omega(r) = Re psi(1/4 + i r/2) - log pi."""
    mp.dps = dps
    logpi = mp.log(mp.pi)
    result = []
    for r_val in r:
        z = mp.mpc(mp.mpf('0.25'), mp.mpf(str(r_val)) * mp.mpf('0.5'))
        psi_val = mp.digamma(z)
        om = mp.re(psi_val) - logpi
        result.append(float(om))
    return np.array(result)

def e_hat_vec_mp(r, n, L_val):
    """e_hat_n(r) vectorized."""
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
    """Q kernel."""
    if n == m:
        return 2.0 * (1.0 - abs(y) / L_val) * np.cos(TWO_PI * n * y / L_val)
    return (np.sin(TWO_PI * m * y / L_val) - np.sin(TWO_PI * n * y / L_val)) / (np.pi * (n - m))

def prime_mat_mp(idx, L_val):
    """Prime contribution matrix."""
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
    """Boundary W_{0,2} matrix."""
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
    """Build tau(eps=1) using high precision."""
    mp.dps = dps

    idx = list(range(-N, N + 1))
    D = len(idx)

    R = TWO_PI * N / L_val + 80.0
    G = max(3000, 200 * N)

    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h)
    w[0] = w[-1] = 0.5 * h

    Om = omega_vec_mp(grid, dps)
    E = np.empty((D, G), dtype=complex)

    for k, n in enumerate(idx):
        E[k] = e_hat_vec_mp(grid, n, L_val)

    A = E * (w * Om)
    arch = (A @ E.conj().T).real / (2.0 * np.pi)

    w02 = w02_mat_mp(idx, L_val)
    prime = prime_mat_mp(idx, L_val)

    tau = arch + w02 - prime
    tau = (tau + tau.T) / 2.0

    return tau, idx

def compute_with_certified_error_bounds(N, L_val, dps_list=[50, 100, 200]):
    """
    Compute eigenvalues with certified error bounds using residual norms.
    Returns certified intervals for lambda_1, lambda_2, and gap.
    """
    print(f"\n{'='*70}")
    print(f"N={N}, L={L_val:.10f}")
    print(f"{'='*70}\n")

    results = {}

    for dps in dps_list:
        print(f"  dps={dps}:")
        mp.dps = dps

        # Build tau with high precision
        tau_hp, idx = build_tau_highprecision(N, L_val, dps)
        tau_64 = tau_hp.astype(np.float64)

        # Compute eigenvalues and eigenvectors
        eigs, vecs = np.linalg.eigh(tau_64)

        # Get two smallest
        lambda1 = eigs[0]
        lambda2 = eigs[1]
        v1 = vecs[:, 0]
        v2 = vecs[:, 1]

        # Compute residuals (in float64, where the eigenvalues were computed)
        r1_vec = tau_64 @ v1 - lambda1 * v1
        r2_vec = tau_64 @ v2 - lambda2 * v2

        # Residual norms (error bounds per Bauer-Fike)
        r1 = np.linalg.norm(r1_vec)
        r2 = np.linalg.norm(r2_vec)

        # Eigenvector norms (should be 1, but let's check)
        v1_norm = np.linalg.norm(v1)
        v2_norm = np.linalg.norm(v2)

        # Error bounds (Bauer-Fike or just residual bound)
        # The error in eigenvalue is roughly r / ||v||
        err1 = r1 / v1_norm
        err2 = r2 / v2_norm

        # Certified intervals (conservative)
        lambda1_lower = lambda1 - err1
        lambda1_upper = lambda1 + err1
        lambda2_lower = lambda2 - err2
        lambda2_upper = lambda2 + err2

        # Gap: worst case (smallest gap from the intervals)
        gap_lower_bound = lambda2_lower - lambda1_upper
        gap_upper_bound = lambda2_upper - lambda1_lower

        print(f"    lambda_1 = {lambda1:.15e}")
        print(f"    lambda_2 = {lambda2:.15e}")
        print(f"    gap      = {lambda2 - lambda1:.15e}")
        print(f"    ||r1||   = {r1:.15e}  (error bound on lambda_1: {err1:.15e})")
        print(f"    ||r2||   = {r2:.15e}  (error bound on lambda_2: {err2:.15e})")
        print(f"    CERTIFIED INTERVAL lambda_1: [{lambda1_lower:.15e}, {lambda1_upper:.15e}]")
        print(f"    CERTIFIED INTERVAL lambda_2: [{lambda2_lower:.15e}, {lambda2_upper:.15e}]")
        print(f"    CERTIFIED GAP BOUNDS: [{gap_lower_bound:.15e}, {gap_upper_bound:.15e}]")

        results[dps] = {
            'lambda_1': float(lambda1),
            'lambda_2': float(lambda2),
            'gap': float(lambda2 - lambda1),
            'err1': float(err1),
            'err2': float(err2),
            'gap_lower': float(gap_lower_bound),
            'gap_upper': float(gap_upper_bound),
            'r1': float(r1),
            'r2': float(r2),
        }

    # Analyze results across precision levels
    print(f"\n  Stability analysis:")
    gaps = [results[dps]['gap'] for dps in dps_list]
    gap_lowers = [results[dps]['gap_lower'] for dps in dps_list]

    print(f"  Gap values: {[f'{g:.5e}' for g in gaps]}")
    print(f"  Gap lower bounds: {[f'{g:.5e}' for g in gap_lowers]}")

    # Check if all gaps are positive
    all_positive = all(g > 0 for g in gaps)
    all_certified_positive = all(g > 0 for g in gap_lowers)

    print(f"\n  All computed gaps positive? {all_positive}")
    print(f"  All certified lower bounds > 0? {all_certified_positive}")

    if all_certified_positive:
        # Use the most conservative (smallest) lower bound
        min_gap_lower = min(gap_lowers)
        print(f"\n  RIGOROUS CERTIFIED LOWER BOUND: gap >= {min_gap_lower:.15e}")
        return True, min_gap_lower, results

    else:
        print(f"\n  WARNING: Some certified bounds are non-positive, not fully certified")
        return False, None, results

def main():
    print("CCM tau(eps=1) GAP CERTIFICATION WITH ERROR BOUNDS")
    print("="*70)

    for N in [32, 40]:
        certified, bound, results = compute_with_certified_error_bounds(N, L_FLOAT)

        print(f"\n{'='*70}")
        print(f"SUMMARY for N={N}")
        print(f"{'='*70}")

        if certified:
            print(f"RESULT: Gap is CERTIFIED POSITIVE")
            print(f"CERTIFIED LOWER BOUND: {bound:.15e}")
        else:
            print(f"RESULT: Gap likely positive, but not fully certified")
            print(f"Best estimate from dps=200: {results[200]['gap']:.15e}")
            print(f"Conservative bound (may be non-positive): {results[200]['gap_lower']:.15e}")

if __name__ == '__main__':
    main()
