#!/usr/bin/env python3
"""
K=6 detailed nullspace analysis — understand structure of nontrivial solutions.
Focus on the high-nullspace-dimension configs (dim >= 8).
"""

import numpy as np
from sympy import symbols, Matrix, Rational, I, pi, sqrt, sin, exp, N as sym_N
from scipy.linalg import svd, null_space

N_freq = 10
L = float(2 * np.pi)

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

def build_hsum_matrix_numeric_k6(z_vals):
    """Build hsum matrix for K=6."""
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

def z_symbolic_to_numeric(z_sym):
    """Convert symbolic z to numeric complex."""
    re_part, im_part = z_sym.as_real_imag()
    return complex(float(sym_N(re_part)), float(sym_N(im_part)))

def analyze_nullspace(M, z_vals, name):
    """Detailed analysis of a single nullspace."""
    print(f"\n{'='*80}")
    print(f"DETAILED ANALYSIS: {name}")
    print(f"{'='*80}")
    print(f"Configuration z-values:")
    for i, z in enumerate(z_vals, 1):
        print(f"  z{i} = {z:.6f}")

    # Compute SVD
    U, s, Vh = svd(M, full_matrices=True)
    rank = np.sum(s > 1e-10)
    dim = 24 - rank

    print(f"\nSVD Analysis:")
    print(f"  Rank (tol=1e-10): {rank}")
    print(f"  Nullspace dimension: {dim}")
    print(f"  Singular values (last 6): {s[-6:]}")

    # Compute nullspace basis
    if dim > 0:
        null_vecs = null_space(M, rcond=1e-10)
        print(f"\nNullspace basis size: {null_vecs.shape}")

        # Check if solution is multi-block or concentrated
        print(f"\nBlock structure analysis (6 blocks, 4 vars each):")
        block_norms = []
        for block_idx in range(6):
            start_col = block_idx * 4
            end_col = start_col + 4
            block_coeffs = null_vecs[start_col:end_col, :]
            block_norm = np.linalg.norm(block_coeffs)
            block_norms.append(block_norm)
            print(f"  Block {block_idx+1}: ‖coeffs‖ = {block_norm:.6e}")

        total_norm = np.linalg.norm(null_vecs)
        rel_block_norms = [bn / total_norm for bn in block_norms]
        print(f"\nRelative block contribution (normalized):")
        for i, rb in enumerate(rel_block_norms, 1):
            print(f"  Block {i}: {rb*100:.1f}%")

        # Count blocks with significant contribution
        sig_threshold = 0.1 * max(rel_block_norms)
        sig_blocks = sum(1 for rb in rel_block_norms if rb > sig_threshold)
        print(f"\nBlocks with significant contribution (>10% of max): {sig_blocks}/6")

        return dim, sig_blocks
    else:
        print("\nTrivial solution only (dim=0)")
        return 0, 0

# High-nullspace configs from K=5 analogs
high_dim_configs = [
    ('[7] Symmetric 4i (K5 Config 7 analog)',
     [-Rational(5,2)+4*I, -Rational(3,2)+4*I, -Rational(1,2)+4*I, Rational(1,2)+4*I, Rational(3,2)+4*I, Rational(5,2)+4*I]),
    ('[9] SmallReal 7i',
     [-Rational(1,50)+7*I, -Rational(3,100)+7*I, -Rational(1,100)+7*I, Rational(1,100)+7*I, Rational(3,100)+7*I, Rational(1,50)+7*I]),
    ('[14] NegReal 6i (K5 Config 14 analog)',
     [-5+6*I, -3+6*I, -1+6*I, 1+6*I, 3+6*I, 5+6*I]),
    ('[16] Wide symmetric 4i',
     [-3+4*I, -2+4*I, -1+4*I, 1+4*I, 2+4*I, 3+4*I]),
]

summary = []
for name, z_vals_sym in high_dim_configs:
    z_vals = [z_symbolic_to_numeric(z) for z in z_vals_sym]
    M = build_hsum_matrix_numeric_k6(z_vals)
    dim, sig_blocks = analyze_nullspace(M, z_vals, name)
    summary.append((name, dim, sig_blocks))

print(f"\n{'='*80}")
print("SUMMARY")
print(f"{'='*80}")
for name, dim, sig_blocks in summary:
    multi_block = "MULTI-BLOCK" if sig_blocks > 1 else "SINGLE-BLOCK"
    print(f"{name:45s} | dim={dim:2d} | {sig_blocks} significant blocks | {multi_block}")
