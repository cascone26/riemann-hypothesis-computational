#!/usr/bin/env python3
"""
Minimal concrete verification: extract K=6 Config 7 solution, print actual numeric values,
and manually evaluate the cubic polynomial at n=10 to debug why Python says 0
when Lean's compiler says the identity fails.
"""

import numpy as np
from scipy.linalg import null_space

N_freq = 10
L_val = 2 * np.pi

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

# Build matrix and extract solution
M = build_hsum_matrix(k6_z_vals, 6, L_val)
null_vecs = null_space(M, rcond=1e-10)
sol_vec = null_vecs[:, 0]

# Extract p, q, r, s
p_vals = [sol_vec[i*4 + 0] for i in range(6)]
q_vals = [sol_vec[i*4 + 1] for i in range(6)]
r_vals = [sol_vec[i*4 + 2] for i in range(6)]
s_vals = [sol_vec[i*4 + 3] for i in range(6)]

print("="*80)
print("K=6 CONFIG 7 EXTRACTED SOLUTION")
print("="*80)
print(f"Solution vector norm: {np.linalg.norm(sol_vec):.6e}")
print(f"hsum constraint residual: {np.linalg.norm(M @ sol_vec):.6e}")
print()

print("Coefficients:")
for i in range(6):
    print(f"Block {i+1}: p={p_vals[i]:.10e}, q={q_vals[i]:.10e}, r={r_vals[i]:.10e}, s={s_vals[i]:.10e}")

print("\n" + "="*80)
print("MANUAL EVALUATION AT n=10 (CENTER FREQUENCY NODE)")
print("="*80)

n_test = 10
z_vals = k6_z_vals
L = L_val

# Compute A_i and cubic_i for each block
sqrt2 = np.sqrt(2)
A_vals = [sqrt2 * np.sin(L * z / 2) for z in z_vals]
A_conj_vals = [np.conj(A) for A in A_vals]

xi = 2 * np.pi * (n_test - N_freq) / L  # = 0 at n=10 since n=N

print(f"n = {n_test}, N = {N_freq}")
print(f"xi = 2π({n_test}-{N_freq})/{L:.4f} = {xi:.10e}")
print()

print("A_i values:")
for i, A in enumerate(A_vals):
    print(f"  A_{i+1} = {A:.10e}")
print()

# Compute cubics
cubics = []
for i in range(6):
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

print("Cubic values:")
for i, c in enumerate(cubics):
    print(f"  cubic_{i+1} = {c:.10e}")
print()

# Compute quartics
quartics = []
for i in range(6):
    z_i = z_vals[i]
    z_i_conj = np.conj(z_i)
    quartic_i = (z_i + xi) * (z_i - xi) * (z_i_conj + xi) * (z_i_conj - xi)
    quartics.append(quartic_i)

print("Quartic values:")
for i, Q in enumerate(quartics):
    print(f"  Q_{i+1} = {Q:.10e}")
print()

# Cross-product sum
total = 0
for i in range(6):
    term = cubics[i]
    for j in range(6):
        if i != j:
            term = term * quartics[j]
    contribution = term
    total += contribution
    print(f"  cubic_{i+1} * (other quartics) = {contribution:.10e}")

print()
print(f"TOTAL CUBIC POLYNOMIAL SUM = {total:.10e}")
print(f"Magnitude: {np.abs(total):.10e}")
print()

if np.abs(total) < 1e-10:
    print("✗ RESULT = ~0 (but Lean says identity FAILS)")
    print("  This contradicts S2699's definitive Lean compilation failure!")
else:
    print(f"✓ RESULT ≠ 0 (magnitude {np.abs(total):.3e})")
    print("  This is consistent with S2699 (Lean compilation failure).")
