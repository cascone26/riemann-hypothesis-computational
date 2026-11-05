#!/usr/bin/env python3
"""
Sanity check: Run the ADAPTED K=6 script in K=5 mode on a known K=5 config.
Should reproduce K=5's original results if adaptation is correct.
"""

import numpy as np
from sympy import Rational, I, pi, N as sym_N
from scipy.linalg import svd

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

def build_hsum_matrix_numeric_generic(z_vals, K):
    """
    Build hsum matrix for generic K blocks.
    z_vals: list of K complex z-values
    K: number of blocks (5 for K=5, 6 for K=6)
    Returns: 21 x (4*K) matrix
    """
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

def estimate_rank(M, tol=1e-10):
    """Estimate rank using SVD."""
    _, s, _ = svd(M)
    return np.sum(s > tol), s

print("=" * 80)
print("SANITY CHECK: Reproduce K=5 results using adapted K=6 script")
print("=" * 80)
print()

# Test configs from K=5's original survey that gave TRIVIAL SOLUTIONS (dim=0)
k5_trivial_configs = [
    ('[1] Simple (1..5) + 5i [K5 ORIGINAL]', [1+5*I, 2+5*I, 3+5*I, 4+5*I, 5+5*I]),
    ('[2] High-Im 10i [K5 ORIGINAL]', [1+10*I, 2+10*I, 3+10*I, 4+10*I, 5+10*I]),
    ('[8] Powers2 5i [K5 ORIGINAL]', [1+5*I, 2+5*I, 4+5*I, 8+5*I, 16+5*I]),
]

# Test configs from K=5's original survey that gave NONTRIVIAL SOLUTIONS (dim>0)
k5_nontrivial_configs = [
    ('[7] Symmetric 4i [K5 ORIGINAL]', [-2+4*I, -1+4*I, 0+4*I, 1+4*I, 2+4*I]),
    ('[9] SmallReal 7i [K5 ORIGINAL]', [-Rational(1,100)+7*I, -Rational(5,100)+7*I, 0+7*I, Rational(5,100)+7*I, Rational(1,100)+7*I]),
]

print("TEST 1: K=5 configs that SHOULD give dim=0 (trivial only)")
print("-" * 80)
for name, z_vals_sym in k5_trivial_configs:
    z_vals = [z_symbolic_to_numeric(z) for z in z_vals_sym]
    M = build_hsum_matrix_numeric_generic(z_vals, K=5)
    rank, _ = estimate_rank(M)
    dim = 20 - rank  # 20 unknowns for K=5
    expected = "dim=0 ✓" if dim == 0 else f"dim={dim} ✗ MISMATCH!"
    print(f"{name:40s} → dim={dim:2d} rank={rank:2d} | {expected}")

print()
print("TEST 2: K=5 configs that SHOULD give dim>0 (nontrivial)")
print("-" * 80)
for name, z_vals_sym in k5_nontrivial_configs:
    z_vals = [z_symbolic_to_numeric(z) for z in z_vals_sym]
    M = build_hsum_matrix_numeric_generic(z_vals, K=5)
    rank, _ = estimate_rank(M)
    dim = 20 - rank  # 20 unknowns for K=5
    expected = "✓ matches" if dim > 0 else "✗ MISMATCH! Should be nontrivial!"
    print(f"{name:40s} → dim={dim:2d} rank={rank:2d} | {expected}")

print()
print("=" * 80)
print("VERIFICATION SUMMARY")
print("=" * 80)
print()

# Now retest K=6 to make sure the generic K-agnostic matrix builder works
print("CROSS-CHECK: Retest a K=6 config with the generic builder")
print("-" * 80)
k6_config = ([Rational(5,2)+4*I, -Rational(3,2)+4*I, -Rational(1,2)+4*I,
              Rational(1,2)+4*I, Rational(3,2)+4*I, Rational(5,2)+4*I],
             "Config 7 Symmetric 4i (K=6, wide version)")

z_vals = [z_symbolic_to_numeric(z) for z in k6_config[0]]
M = build_hsum_matrix_numeric_generic(z_vals, K=6)
rank, _ = estimate_rank(M)
dim = 24 - rank
print(f"K=6 Config 7 Symmetric 4i (wide) → dim={dim:2d} rank={rank:2d}")
print(f"  (Expected: dim>0, ideally dim≥8 for good multi-block solution)")

print()
print("CONCLUSION:")
print("  If K=5 trivial configs still give dim=0 and K=5 nontrivial still give dim>0,")
print("  then the adapted script is CORRECT and K=6's generic solvability is REAL.")
print("  If K=5 trivial configs now give dim>0, then there is a BUG in the adaptation.")
