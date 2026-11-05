#!/usr/bin/env python3
"""
Verify that the extracted nullspace solutions actually satisfy the hsum constraint.
This is critical: if they don't, the polynomial evaluation result is meaningless.
"""

import numpy as np
from scipy.linalg import null_space

N_freq = 10
L_val = 2 * np.pi

# K=5 Config 7
k5_z_vals = [-2+4j, -1+4j, 0+4j, 1+4j, 2+4j]

# K=6 Config 7
k6_z_vals = [-2.5+4j, -1.5+4j, -0.5+4j, 0.5+4j, 1.5+4j, 2.5+4j]

def v_plus_lc_numeric(n, L, z):
    """Evaluate vPlusLC numerically."""
    freq_shift = 2 * np.pi * (n - N_freq) / L
    denom = z + freq_shift
    num = np.sqrt(2) * np.sin(L * z / 2)
    return num / denom

def v_minus_lc_numeric(n, L, z):
    """Evaluate vMinusLC numerically."""
    freq_shift = 2 * np.pi * (n - N_freq) / L
    denom = z - freq_shift
    num = np.sqrt(2) * np.sin(L * z / 2)
    return num / denom

def build_hsum_matrix(z_vals, K, L):
    """Build hsum constraint matrix for K blocks."""
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

def verify_nullspace(K, config_name, z_vals):
    """Verify that the extracted nullspace solution satisfies M*x=0."""
    print(f"\n{'='*80}")
    print(f"VERIFICATION: K={K} {config_name}")
    print(f"{'='*80}")

    M = build_hsum_matrix(z_vals, K, L_val)
    null_vecs = null_space(M, rcond=1e-10)

    print(f"Matrix shape: {M.shape} ({2*N_freq+1} equations, {4*K} unknowns)")
    print(f"Nullspace dimension: {null_vecs.shape[1]}")

    if null_vecs.shape[1] == 0:
        print(f"ERROR: No nullspace found!")
        return False

    # Test the first basis vector
    x = null_vecs[:, 0]
    residual = M @ x
    residual_norm = np.linalg.norm(residual)

    print(f"\nNull-space vector norm: {np.linalg.norm(x):.6e}")
    print(f"Residual M*x norm: {residual_norm:.6e}")
    print(f"Max residual component: {np.max(np.abs(residual)):.6e}")

    # Check constraint per frequency node
    print(f"\nVerification across all 21 frequency nodes:")
    print(f"{'n':3s} | {'|residual(n)|':20s} | Status")
    print("-" * 50)

    all_ok = True
    for n in range(2*N_freq+1):
        row_start = n * 4 * K
        row_end = row_start + 4 * K
        node_residual = np.linalg.norm(residual[row_start:row_end])
        status = "✓" if node_residual < 1e-10 else "✗"
        print(f"{n:3d} | {node_residual:20.6e} | {status}")
        if node_residual >= 1e-10:
            all_ok = False

    if all_ok:
        print(f"\n✓ All frequency nodes satisfy hsum constraint (residual < 1e-10)")
        return True
    else:
        print(f"\n✗ Some frequency nodes DO NOT satisfy hsum constraint")
        return False

if __name__ == "__main__":
    print("\n" + "="*80)
    print("HSUM CONSTRAINT VERIFICATION")
    print("="*80)

    k5_ok = verify_nullspace(5, "Config 7 (Symmetric -2..2+4i)", k5_z_vals)
    k6_ok = verify_nullspace(6, "Config 7 (Symmetric -2.5..2.5+4i)", k6_z_vals)

    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print(f"K=5 constraint satisfied: {k5_ok}")
    print(f"K=6 constraint satisfied: {k6_ok}")
    print(f"{'='*80}\n")
