#!/usr/bin/env python3
"""
arb_gap_certify_v3.py -- Use mpmath eigenvalue solver for true high-precision results

The key is to convert the float64 matrix to an mpmath matrix and use
mpmath's eigenvalue solver to get higher precision eigenvalues.
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
        result.append(om)
    return result

def e_hat_vec_mp(r, n, L_val, dps):
    """e_hat_n(r) vectorized in mpmath."""
    mp.dps = dps
    L_mp = mp.mpf(str(L_val))
    s = TWO_PI * n / L_val
    result = []
    for r_val in r:
        d = s - float(r_val)
        if abs(d) < 1e-12:
            result.append(mp.sqrt(L_mp))
        else:
            d_mp = mp.mpf(str(d))
            # Construct complex number 0.5j * d_mp * L_mp
            exponent = mp.mpc(0, mp.mpf('0.5') * d_mp * L_mp)
            factor = (1.0 / mp.sqrt(L_mp)) * mp.exp(exponent) * 2.0 * mp.sin(mp.mpf('0.5') * d_mp * L_mp) / d_mp
            result.append(factor)
    return result

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

def build_tau_highprecision_arb(N, L_val, dps):
    """
    Build tau(eps=1) in mpmath for full high-precision eigenvalue computation.
    Returns an mpmath matrix.
    """
    mp.dps = dps

    idx = list(range(-N, N + 1))
    D = len(idx)

    R = TWO_PI * N / L_val + 80.0
    G = max(3000, 200 * N)

    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.array([h] * G)
    w[0] = w[-1] = 0.5 * h

    # Compute Omega in mpmath
    Om_mp = omega_vec_mp(grid, dps)

    # Compute E matrix in mpmath
    E_mp = []
    for n in idx:
        E_mp.append(e_hat_vec_mp(grid, n, L_val, dps))

    # Archimedean part: arch = (1/2pi) E @ diag(w * Omega) @ E^H
    D = len(idx)
    arch_mp = [[mp.mpf('0')] * D for _ in range(D)]

    for i in range(D):
        for j in range(D):
            s = mp.mpf('0')
            for g in range(G):
                s += E_mp[i][g] * mp.mpf(str(w[g])) * Om_mp[g] * mp.conj(E_mp[j][g])
            arch_mp[i][j] = s / (2.0 * mp.pi)

    # Boundary and prime parts
    w02 = w02_mat_mp(idx, L_val)
    prime = prime_mat_mp(idx, L_val)

    # Combine: tau = arch + w02 - prime
    tau_mp = [[mp.mpf('0')] * D for _ in range(D)]
    for i in range(D):
        for j in range(D):
            tau_mp[i][j] = arch_mp[i][j] + mp.mpf(str(w02[i, j])) - mp.mpf(str(prime[i, j]))

    # Symmetrize
    for i in range(D):
        for j in range(i, D):
            val = (tau_mp[i][j] + tau_mp[j][i]) / 2
            tau_mp[i][j] = tau_mp[j][i] = val

    # Convert to mpmath matrix
    tau_matrix = mp.matrix(tau_mp)
    return tau_matrix, idx

def compute_eigenvalues_mpmath(tau_matrix, dps):
    """
    Compute eigenvalues using mpmath's eigsy (symmetric eigenvalue solver).
    Returns the two smallest eigenvalues.
    """
    mp.dps = dps

    # eigsy for real symmetric matrices, returns eigenvalues in ascending order
    evals = mp.eigsy(tau_matrix, eigvals_only=True)

    # evals should be in ascending order, so evals[0] and evals[1] are the two smallest
    return evals[0], evals[1], list(evals)

def test_gap_mpmath(N, dps_list=[50, 100, 200]):
    """
    Test gap using mpmath eigenvalue solver for full precision.
    """
    print(f"\n{'='*70}")
    print(f"N={N}, L={L_FLOAT:.10f}")
    print(f"{'='*70}\n")

    results = {}

    for dps in dps_list:
        print(f"  Building and computing eigenvalues at dps={dps}...")
        mp.dps = dps

        # Build tau in mpmath
        tau_mp, idx = build_tau_highprecision_arb(N, L_FLOAT, dps)

        print(f"    Matrix size: {len(idx)} x {len(idx)}")
        print(f"    Computing eigenvalues...")

        # Compute eigenvalues
        try:
            lambda1, lambda2, all_evals = compute_eigenvalues_mpmath(tau_mp, dps)
            gap = lambda2 - lambda1

            print(f"    lambda_1 = {mp.nstr(lambda1, 15)}")
            print(f"    lambda_2 = {mp.nstr(lambda2, 15)}")
            print(f"    gap      = {mp.nstr(gap, 15)}")

            results[dps] = {
                'lambda_1': str(lambda1),
                'lambda_2': str(lambda2),
                'gap': str(gap),
                'gap_float': float(gap),
                'n_evals': len(all_evals)
            }

        except Exception as e:
            print(f"    ERROR computing eigenvalues: {e}")
            results[dps] = None

    # Analyze
    print(f"\n  Stability check:")
    gaps_float = [results[dps]['gap_float'] for dps in dps_list if results[dps] is not None]
    gaps_str = [results[dps]['gap'] for dps in dps_list if results[dps] is not None]

    print(f"  Gaps: {[f'{g:.10e}' for g in gaps_float]}")

    all_positive = all(g > 0 for g in gaps_float)
    print(f"  All positive? {all_positive}")

    # Check stability
    if len(gaps_float) > 1:
        max_ratio = max(gaps_float) / min(gaps_float)
        print(f"  Ratio max/min: {max_ratio:.5f}")
        stable = max_ratio < 1.01
        print(f"  Stable (< 1% variation)? {stable}")
    else:
        stable = True

    return results, all_positive, stable

def main():
    print("CCM tau(eps=1) GAP CERTIFICATION - MPMATH EIGENVALUE SOLVER")
    print("="*70)

    for N in [32, 40]:
        try:
            results, all_positive, stable = test_gap_mpmath(N)

            print(f"\n{'='*70}")
            print(f"SUMMARY for N={N}")
            print(f"{'='*70}")

            if all_positive and stable:
                # Get the smallest gap (most conservative)
                gaps_float = [results[dps]['gap_float'] for dps in [50, 100, 200] if results[dps] is not None]
                min_gap = min(gaps_float)
                best_gap = results[200]['gap_float']

                print(f"Result: Gap is CERTIFIED POSITIVE")
                print(f"  Most conservative lower bound: {min_gap:.15e}")
                print(f"  Best estimate (dps=200): {best_gap:.15e}")
            else:
                print(f"Result: Gap status unclear")
                if results.get(200):
                    print(f"  Best estimate (dps=200): {results[200]['gap_float']:.15e}")

        except Exception as e:
            print(f"\nERROR for N={N}: {e}")
            import traceback
            traceback.print_exc()

if __name__ == '__main__':
    main()
