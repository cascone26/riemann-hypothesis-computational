#!/usr/bin/env python3
"""
S2742 Verification: High-precision numeric evaluation of K=5 Config 7 cubic polynomial.
Uses mpmath (dps=50+), NOT sympy.simplify().
"""

import mpmath as mp
import numpy as np
from scipy.linalg import null_space

# High precision setup
mp.mp.dps = 60  # 60 decimal places

N_freq = 10
L_val_float = 2 * np.pi

# S2742 K=5 Config 7: Symmetric -2..2 + 4i
k5_z_vals = [-2+4j, -1+4j, 0+4j, 1+4j, 2+4j]

def v_plus_lc_numeric(n, L, z):
    """Numeric evaluation of vPlusLC."""
    freq_shift = 2 * np.pi * (n - N_freq) / L
    denom = z + freq_shift
    num = np.sqrt(2) * np.sin(L * z / 2)
    return num / denom

def v_minus_lc_numeric(n, L, z):
    """Numeric evaluation of vMinusLC."""
    freq_shift = 2 * np.pi * (n - N_freq) / L
    denom = z - freq_shift
    num = np.sqrt(2) * np.sin(L * z / 2)
    return num / denom

def build_hsum_matrix(z_vals, K, L):
    """Build hsum matrix."""
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

# Extract K=5 Config 7 solution
print("="*80)
print("S2742 VERIFICATION: K=5 Config 7 High-Precision Evaluation")
print("="*80)
print(f"Configuration: z = {{-2+4i, -1+4i, 0+4i, 1+4i, 2+4i}}")
print(f"Precision: mpmath dps={mp.mp.dps}")
print()

M = build_hsum_matrix(k5_z_vals, 5, L_val_float)
null_vecs = null_space(M, rcond=1e-10)
sol_vec = null_vecs[:, 0]

# Verify hsum constraint
hsum_residual = np.linalg.norm(M @ sol_vec)
print(f"hsum constraint residual: {hsum_residual:.6e}")
print(f"Nullspace dimension: {null_vecs.shape[1]}")
print()

# Extract coefficients
p_vals = [sol_vec[i*4 + 0] for i in range(5)]
q_vals = [sol_vec[i*4 + 1] for i in range(5)]
r_vals = [sol_vec[i*4 + 2] for i in range(5)]
s_vals = [sol_vec[i*4 + 3] for i in range(5)]

# High-precision cubic evaluation at each frequency node
print("High-precision cubic evaluation (mpmath dps=60):")
print(f"{'n':3s} | {'Real part':25s} | {'Imag part':25s} | {'Magnitude':20s}")
print("-" * 100)

results = []
for n_test in range(2 * N_freq + 1):
    # Convert to mpmath complex
    L_mp = mp.mpf('2') * mp.pi
    xi_mp = mp.mpf('2') * mp.pi * (mp.mpf(str(n_test)) - mp.mpf(str(N_freq))) / L_mp

    z_mp = [mp.mpc(mp.mpf(str(z.real)), mp.mpf(str(z.imag))) for z in k5_z_vals]
    p_mp = [mp.mpc(mp.mpf(str(p.real)), mp.mpf(str(p.imag))) for p in p_vals]
    q_mp = [mp.mpc(mp.mpf(str(q.real)), mp.mpf(str(q.imag))) for q in q_vals]
    r_mp = [mp.mpc(mp.mpf(str(r.real)), mp.mpf(str(r.imag))) for r in r_vals]
    s_mp = [mp.mpc(mp.mpf(str(s.real)), mp.mpf(str(s.imag))) for s in s_vals]

    # Compute A_i (using mpmath sin)
    sqrt2_mp = mp.sqrt(mp.mpf('2'))
    A_vals_mp = [sqrt2_mp * mp.sin(L_mp * z / mp.mpf('2')) for z in z_mp]
    A_conj_vals_mp = [mp.conj(A) for A in A_vals_mp]

    # Compute cubics
    cubics_mp = []
    for i in range(5):
        z_i = z_mp[i]
        z_i_conj = mp.conj(z_i)
        A_i = A_vals_mp[i]
        A_i_conj = A_conj_vals_mp[i]
        p_i, q_i, r_i, s_i = p_mp[i], q_mp[i], r_mp[i], s_mp[i]

        cubic_i = (
            (p_i + q_i) * A_i * (z_i * (z_i_conj ** 2))
            + (r_i + s_i) * A_i_conj * ((z_i ** 2) * z_i_conj)
            + ((q_i - p_i) * A_i * (z_i_conj ** 2) + (s_i - r_i) * A_i_conj * (z_i ** 2)) * xi_mp
            + (-(p_i + q_i) * A_i * z_i - (r_i + s_i) * A_i_conj * z_i_conj) * (xi_mp ** 2)
            + ((p_i - q_i) * A_i + (r_i - s_i) * A_i_conj) * (xi_mp ** 3)
        )
        cubics_mp.append(cubic_i)

    # Compute quartics
    quartics_mp = []
    for i in range(5):
        z_i = z_mp[i]
        z_i_conj = mp.conj(z_i)
        quartic_i = (z_i + xi_mp) * (z_i - xi_mp) * (z_i_conj + xi_mp) * (z_i_conj - xi_mp)
        quartics_mp.append(quartic_i)

    # Cross-product sum
    total_mp = mp.mpc('0')
    for i in range(5):
        term = cubics_mp[i]
        for j in range(5):
            if i != j:
                term = term * quartics_mp[j]
        total_mp += term

    # Extract real and imag parts
    real_part = mp.re(total_mp)
    imag_part = mp.im(total_mp)
    magnitude = mp.fabs(total_mp)

    results.append((n_test, total_mp, magnitude))

    # Format output (show full precision for critical nodes)
    if n_test in [0, 5, 10, 15, 20]:
        print(f"{n_test:3d} | {str(real_part):25s} | {str(imag_part):25s} | {str(magnitude):20s}")
    else:
        # Abbreviated for non-critical nodes
        print(f"{n_test:3d} | {mp.nstr(real_part, 15):25s} | {mp.nstr(imag_part, 15):25s} | {mp.nstr(magnitude, 15):20s}")

print()
print("="*80)
print("CRITICAL NODES (Full Precision)")
print("="*80)

for n_test in [0, 5, 10, 15, 20]:
    _, total_mp, _ = results[n_test]
    print(f"\nn={n_test}:")
    print(f"  Real:  {mp.nstr(mp.re(total_mp), 50)}")
    print(f"  Imag:  {mp.nstr(mp.im(total_mp), 50)}")
    print(f"  |val|: {mp.nstr(mp.fabs(total_mp), 50)}")

print()
print("="*80)
print("SUMMARY")
print("="*80)

nonzero_count = sum(1 for _, _, mag in results if mag > mp.mpf('1e-40'))
print(f"Nodes with magnitude > 1e-40: {nonzero_count}/21")

if nonzero_count == 0:
    print("\n✓ K=5 Config 7 CUBIC = 0 (confirmed high-precision)")
    print("  S2742 is CORRECT")
else:
    min_mag = min(mag for _, _, mag in results)
    max_mag = max(mag for _, _, mag in results)
    print(f"\n✗ K=5 Config 7 CUBIC ≠ 0")
    print(f"  Magnitude range: [{mp.nstr(min_mag, 15)} .. {mp.nstr(max_mag, 15)}]")
    print(f"  S2742 is INCORRECT (sympy.simplify() false positive)")
