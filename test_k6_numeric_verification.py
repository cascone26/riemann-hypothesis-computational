#!/usr/bin/env python3
"""
K=6 polynomial evaluation using NUMERIC (float) arithmetic only, not symbolic simplification.
This avoids the sympy.simplify() error that claimed all results = 0.
"""

import numpy as np
from scipy.linalg import null_space

N_freq = 10
L_val = 2 * np.pi

# K=5 Config 7 for sanity check
k5_z_vals = [-2+4j, -1+4j, 0+4j, 1+4j, 2+4j]

# K=6 Config 7 for real test
k6_z_vals = [-2.5+4j, -1.5+4j, -0.5+4j, 0.5+4j, 1.5+4j, 2.5+4j]

def v_plus_lc_numeric(n, L, z):
    freq_shift = 2 * np.pi * (n - N_freq) / L
    denom = z + freq_shift
    num = np.sqrt(2) * np.sin(L * z / 2)
    return num / denom

def v_minus_lc_numeric(n, L, z):
    freq_shift = 2 * np.pi * (n - N_freq) / L
    denom = z - freq_shift
    num = np.sqrt(2) * np.sin(L * z / 2)
    return num / denom

def build_hsum_matrix(z_vals, K, L):
    num_nodes = 2 * N_freq + 1
    M = []
    for n_val in range(num_nodes):
        row = []
        for z_val in z_vals:
            z_conj = np.conj(z_val)
            row.extend([v_plus_lc_numeric(n_val, L, z_val),
                       v_minus_lc_numeric(n_val, L, z_val),
                       v_plus_lc_numeric(n_val, L, z_conj),
                       v_minus_lc_numeric(n_val, L, z_conj)])
        M.append(row)
    return np.array(M, dtype=complex)

def evaluate_polynomial_numeric(K, n_val, N, L, z_vals, p_vals, q_vals, r_vals, s_vals):
    """Evaluate polynomial using pure numeric (float) arithmetic."""
    sqrt2 = np.sqrt(2)
    xi = 2 * np.pi * (n_val - N) / L

    # A_i values
    A_vals = [sqrt2 * np.sin(L * z / 2) for z in z_vals]
    A_conj_vals = [np.conj(A) for A in A_vals]

    # Cubics
    cubics = []
    for i in range(K):
        z_i = z_vals[i]
        z_i_conj = np.conj(z_i)
        A_i = A_vals[i]
        A_i_conj = A_conj_vals[i]
        p_i, q_i, r_i, s_i = p_vals[i], q_vals[i], r_vals[i], s_vals[i]

        cubic_i = (
            (p_i + q_i) * A_i * (z_i * (z_i_conj ** 2))
            + (r_i + s_i) * A_i_conj * ((z_i ** 2) * z_i_conj)
            + ((q_i - p_i) * A_i * (z_i_conj ** 2) + (s_i - r_i) * A_i_conj * (z_i ** 2)) * xi
            + (-(p_i + q_i) * A_i * z_i - (r_i + s_i) * A_i_conj * z_i_conj) * (xi ** 2)
            + ((p_i - q_i) * A_i + (r_i - s_i) * A_i_conj) * (xi ** 3)
        )
        cubics.append(cubic_i)

    # Quartics
    quartics = []
    for i in range(K):
        z_i = z_vals[i]
        z_i_conj = np.conj(z_i)
        quartic_i = (z_i + xi) * (z_i - xi) * (z_i_conj + xi) * (z_i_conj - xi)
        quartics.append(quartic_i)

    # Cross-product sum
    total = 0
    for i in range(K):
        term = cubics[i]
        for j in range(K):
            if i != j:
                term = term * quartics[j]
        total += term

    return total

def test_config(K, name, z_vals):
    """Test one K configuration with numeric evaluation."""
    print(f"\n{'='*80}")
    print(f"K={K} {name}")
    print(f"{'='*80}")

    # Extract nullspace solution
    M = build_hsum_matrix(z_vals, K, L_val)
    null_vecs = null_space(M, rcond=1e-10)
    sol_vec = null_vecs[:, 0]

    p_vals = [sol_vec[i*4 + 0] for i in range(K)]
    q_vals = [sol_vec[i*4 + 1] for i in range(K)]
    r_vals = [sol_vec[i*4 + 2] for i in range(K)]
    s_vals = [sol_vec[i*4 + 3] for i in range(K)]

    print(f"Nullspace dimension: {null_vecs.shape[1]}")
    print(f"hsum constraint residual: {np.linalg.norm(M @ sol_vec):.3e}")
    print()

    # Evaluate polynomial at all 21 frequency nodes (NUMERIC ONLY)
    print(f"{'n':3s} | {'Magnitude':20s} | Status")
    print("-" * 50)

    results = []
    for n_test in range(2 * N_freq + 1):
        poly = evaluate_polynomial_numeric(K, n_test, N_freq, L_val, z_vals, p_vals, q_vals, r_vals, s_vals)
        mag = np.abs(poly)
        is_zero = mag < 1e-8  # Numeric zero threshold
        status = "✓ ZERO" if is_zero else f"✗ NONZERO"
        results.append((n_test, mag, is_zero))
        print(f"{n_test:3d} | {mag:20.6e} | {status}")

    print()
    zero_count = sum(1 for _, _, is_z in results if is_z)
    total_count = len(results)
    print(f"Result: {zero_count}/{total_count} zeros")

    if zero_count == total_count:
        print(f"✓ K={K} POLYNOMIAL = 0 for all nodes")
        return "PASS"
    else:
        nonzero_count = total_count - zero_count
        avg_nonzero = np.mean([mag for _, mag, is_z in results if not is_z])
        print(f"✗ K={K} POLYNOMIAL ≠ 0 for {nonzero_count} nodes (avg magnitude {avg_nonzero:.3e})")
        return "FAIL"

if __name__ == "__main__":
    print("\n" + "="*80)
    print("K=5 AND K=6 NUMERIC POLYNOMIAL EVALUATION (NO SYMBOLIC SIMPLIFICATION)")
    print("="*80)

    result_k5 = test_config(5, "Config 7 (Symmetric -2..2+4i, sanity check)", k5_z_vals)
    result_k6 = test_config(6, "Config 7 (Symmetric -2.5..2.5+4i, real test)", k6_z_vals)

    print(f"\n{'='*80}")
    print("FINAL RESULTS (NUMERIC, NOT SYMBOLIC)")
    print(f"{'='*80}")
    print(f"K=5 Config 7: {result_k5}")
    print(f"K=6 Config 7: {result_k6}")
    print(f"{'='*80}\n")
