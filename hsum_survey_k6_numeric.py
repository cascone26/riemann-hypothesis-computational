#!/usr/bin/env python3
"""
K=6 hsum constraint nullspace survey — NUMERIC (fast)
Evaluates symbolic expressions to high-precision floats, computes rank numerically.
Uses singular value decomposition to estimate rank. N=10 (full constraint).
"""

import numpy as np
from sympy import symbols, Matrix, Rational, I, pi, sqrt, sin, conjugate, exp, N as sym_N
from scipy.linalg import svd

N_freq = 10  # Full constraint
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
    """Build hsum matrix numerically for K=6 (6 blocks, 24 unknowns)."""
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

def estimate_rank(M, tol=1e-10):
    """Estimate rank using SVD."""
    _, s, _ = svd(M)
    return np.sum(s > tol), s

def z_symbolic_to_numeric(z_sym):
    """Convert symbolic z to numeric complex."""
    re_part, im_part = z_sym.as_real_imag()
    return complex(float(sym_N(re_part)), float(sym_N(im_part)))

# K=6 configurations: try analogs of K=5 successes plus new ones
configs = [
    ('[1] Simple (1..6) + 5i', [1+5*I, 2+5*I, 3+5*I, 4+5*I, 5+5*I, 6+5*I]),
    ('[2] High-Im 10i', [1+10*I, 2+10*I, 3+10*I, 4+10*I, 5+10*I, 6+10*I]),
    ('[3] Clustered 8i', [Rational(1)+8*I, Rational(11,10)+8*I, Rational(12,10)+8*I, Rational(13,10)+8*I, Rational(14,10)+8*I, Rational(15,10)+8*I]),
    ('[4] Spread 3i', [Rational(1,2)+3*I, 2+3*I, Rational(7,2)+3*I, 5+3*I, 7+3*I, 9+3*I]),
    ('[5] π-mult 3i', [pi+3*I, 2*pi+3*I, 3*pi+3*I, 4*pi+3*I, 5*pi+3*I, 6*pi+3*I]),
    ('[6] VeryHigh 20i', [1+20*I, 2+20*I, 3+20*I, 4+20*I, 5+20*I, 6+20*I]),
    # K=5 Config 7 analog: Symmetric real parts spanning a range, constant Im
    ('[7] Symmetric 4i (K5 Config 7 analog)', [-Rational(5,2)+4*I, -Rational(3,2)+4*I, -Rational(1,2)+4*I, Rational(1,2)+4*I, Rational(3,2)+4*I, Rational(5,2)+4*I]),
    ('[8] Powers2 5i', [1+5*I, 2+5*I, 4+5*I, 8+5*I, 16+5*I, 32+5*I]),
    ('[9] SmallReal 7i', [-Rational(1,50)+7*I, -Rational(3,100)+7*I, -Rational(1,100)+7*I, Rational(1,100)+7*I, Rational(3,100)+7*I, Rational(1,50)+7*I]),
    ('[10] VaryIm', [1+2*I, 2+3*I, 3+4*I, 4+5*I, 5+6*I, 6+7*I]),
    ('[11] e-mult 5i', [exp(1)+5*I, 2*exp(1)+5*I, 3*exp(1)+5*I, 4*exp(1)+5*I, 5*exp(1)+5*I, 6*exp(1)+5*I]),
    ('[12] Golden φ^k', [((1+sqrt(5))/2)**k + 5*I for k in range(6)]),
    ('[13] Irregular', [Rational(3,10)+Rational(21,10)*I, Rational(17,10)+Rational(23,10)*I,
                       Rational(9,10)+Rational(26,10)*I, Rational(12,5)+Rational(11,5)*I,
                       Rational(11,10)+Rational(14,5)*I, Rational(2,1)+Rational(8,5)*I]),
    ('[14] NegReal 6i (K5 Config 14 analog)', [-5+6*I, -3+6*I, -1+6*I, 1+6*I, 3+6*I, 5+6*I]),
    ('[15] MixedRat 5i', [Rational(1,3)+5*I, Rational(2,3)+5*I, Rational(4,3)+5*I, Rational(5,3)+5*I, 2+5*I, Rational(7,3)+5*I]),
    # Additional configurations for K=6 specifically
    ('[16] Wide symmetric 4i', [-3+4*I, -2+4*I, -1+4*I, 1+4*I, 2+4*I, 3+4*I]),
    ('[17] Small-real cluster 7i', [-Rational(2,100)+7*I, -Rational(1,100)+7*I, 0+7*I, Rational(1,100)+7*I, Rational(2,100)+7*I, Rational(3,100)+7*I]),
]

print("K=6 hsum Nullspace Dimension Survey (NUMERIC)")
print("=" * 80)
print(f"N = {N_freq} ({2*N_freq+1} equations, FULL constraint)")
print(f"K = 6 blocks (24 unknowns: p_i, q_i, r_i, s_i for i=1..6)")
print(f"Matrix shape: 21 × 24")
print(f"Method: SVD on numerically evaluated matrices (tol=1e-10)")
print("=" * 80)
print()

results = []
for name, z_vals_sym in configs:
    print(f"{name:45s}", end=' ', flush=True)
    try:
        # Convert to numeric
        z_vals = [z_symbolic_to_numeric(z) for z in z_vals_sym]

        # Build and analyze matrix
        M = build_hsum_matrix_numeric_k6(z_vals)
        rank, s_values = estimate_rank(M)
        dim = 24 - rank  # 24 unknowns for K=6

        print(f"rank={rank:2d} dim={dim:2d}")
        results.append((name, dim, rank))
    except Exception as e:
        print(f"ERROR: {str(e)[:40]}")

print()
print("=" * 80)
dims = [r[1] for r in results]
if dims:
    print(f"Total: {len(results)} configs tested")
    print(f"Nullspace dims observed: {sorted(set(dims))}")
    print(f"Configs with dim=0 (trivial only): {sum(1 for d in dims if d == 0)}")
    print(f"Configs with dim>0 (nontrivial): {sum(1 for d in dims if d > 0)}")
    print()
    if all(d == 0 for d in dims):
        print("*** PATTERN: ALL CONFIGS HAVE TRIVIAL SOLUTIONS ONLY (dim=0) ***")
    elif all(d > 0 for d in dims):
        print("*** PATTERN: ALL CONFIGS HAVE NONTRIVIAL SOLUTIONS (dim>0) ***")
    else:
        zero_count = sum(1 for d in dims if d == 0)
        pos_count = sum(1 for d in dims if d > 0)
        print(f"*** PATTERN: MIXED — {zero_count} zero, {pos_count} positive ***")

print()
print("Detailed Results:")
print()
for name, dim, rank in results:
    status = "Nontrivial" if dim > 0 else "Trivial only"
    print(f"{name:45s} | dim={dim:2d} rank={rank:2d} | {status}")
