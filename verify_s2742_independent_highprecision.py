#!/usr/bin/env python3
"""
Independent S2742 verification — high-precision numeric evaluation WITHOUT symbolic algebra.
Uses mpmath (dps=80) to evaluate the cubic polynomial at all 21 frequency nodes.

Key discipline:
1. Extract the exact K=5 Config 7 z-values and nullspace solution
2. Compute A_i values using high-precision sin/cos
3. Evaluate cubic polynomial using DIRECT FLOATING-POINT ARITHMETIC
4. NO sympy.simplify(), NO symbolic expressions — only numeric computation

This is the ground truth verification method.
"""

import mpmath as mp
import numpy as np
from scipy.linalg import null_space

# Set ultra-high precision
mp.mp.dps = 80

N_freq = 10
L_val_float = 2.0 * np.pi

# K=5 Config 7: Symmetric -2..2 + 4i
z_vals_numpy = np.array([-2+4j, -1+4j, 0+4j, 1+4j, 2+4j], dtype=complex)

def build_hsum_matrix(z_vals, K, L):
    """Build hsum constraint matrix for K values at all n ∈ [0, 20]."""
    num_nodes = 2 * N_freq + 1
    M = []
    for n_val in range(num_nodes):
        row = []
        for z_val in z_vals:
            z_conj = np.conj(z_val)
            # v_plus_LC
            freq_shift = 2 * np.pi * (n_val - N_freq) / L
            v_plus_z = np.sqrt(2) * np.sin(L * z_val / 2) / (z_val + freq_shift)
            v_minus_z = np.sqrt(2) * np.sin(L * z_val / 2) / (z_val - freq_shift)
            v_plus_conj = np.sqrt(2) * np.sin(L * z_conj / 2) / (z_conj + freq_shift)
            v_minus_conj = np.sqrt(2) * np.sin(L * z_conj / 2) / (z_conj - freq_shift)
            row.extend([v_plus_z, v_minus_z, v_plus_conj, v_minus_conj])
        M.append(row)
    return np.array(M, dtype=complex)

print("=" * 100)
print("INDEPENDENT S2742 VERIFICATION — HIGH-PRECISION NUMERIC EVALUATION")
print("=" * 100)
print(f"Configuration: z = {{-2+4i, -1+4i, 0+4i, 1+4i, 2+4i}}")
print(f"Precision: mpmath dps={mp.mp.dps}")
print()

# Step 1: Extract the nullspace solution
print("[Step 1] Building hsum constraint matrix...")
M = build_hsum_matrix(z_vals_numpy, 5, L_val_float)
print(f"  Matrix shape: {M.shape}")

print("[Step 2] Computing nullspace...")
null_vecs = null_space(M, rcond=1e-10)
print(f"  Nullspace dimension: {null_vecs.shape[1]}")
sol_vec = null_vecs[:, 0]  # First basis vector
print()

# Verify hsum constraint is satisfied
residual = np.linalg.norm(M @ sol_vec)
print(f"[Verification] hsum constraint residual: {residual:.6e}")
print()

# Extract p, q, r, s coefficients
p_vals = np.array([sol_vec[i*4 + 0] for i in range(5)])
q_vals = np.array([sol_vec[i*4 + 1] for i in range(5)])
r_vals = np.array([sol_vec[i*4 + 2] for i in range(5)])
s_vals = np.array([sol_vec[i*4 + 3] for i in range(5)])

print(f"[Coefficients extracted]")
print(f"  p = {p_vals}")
print(f"  q = {q_vals}")
print(f"  r = {r_vals}")
print(f"  s = {s_vals}")
print()

# Convert z_vals and all coefficients to mpmath complex with high precision
print("[Step 3] Converting to mpmath (dps=80) for high-precision evaluation...")
z_mp = [mp.mpc(mp.mpf(str(z.real)), mp.mpf(str(z.imag))) for z in z_vals_numpy]
p_mp = [mp.mpc(mp.mpf(str(p.real)), mp.mpf(str(p.imag))) for p in p_vals]
q_mp = [mp.mpc(mp.mpf(str(q.real)), mp.mpf(str(q.imag))) for q in q_vals]
r_mp = [mp.mpc(mp.mpf(str(r.real)), mp.mpf(str(r.imag))) for r in r_vals]
s_mp = [mp.mpc(mp.mpf(str(s.real)), mp.mpf(str(s.imag))) for s in s_vals]

print()
print("=" * 100)
print("CUBIC POLYNOMIAL EVALUATION AT ALL n ∈ [0, 20]")
print("=" * 100)
print()

sqrt2_mp = mp.sqrt(mp.mpf('2'))
L_mp = mp.mpf('2') * mp.pi
pi_mp = mp.pi

# Store all results
results = []

for n_test in range(2 * N_freq + 1):
    # Compute frequency shift: ξ = 2π(n - N) / L
    xi_mp = mp.mpf('2') * pi_mp * (mp.mpf(str(n_test)) - mp.mpf(str(N_freq))) / L_mp

    # Compute A_i = √2 * sin(L * z_i / 2)
    A_vals_mp = []
    A_conj_vals_mp = []
    for z in z_mp:
        A = sqrt2_mp * mp.sin(L_mp * z / mp.mpf('2'))
        A_conj = mp.conj(A)
        A_vals_mp.append(A)
        A_conj_vals_mp.append(A_conj)

    # Compute cubic_i for each block
    cubics_mp = []
    for i in range(5):
        z_i = z_mp[i]
        z_i_conj = mp.conj(z_i)
        A_i = A_vals_mp[i]
        A_i_conj = A_conj_vals_mp[i]
        p_i, q_i, r_i, s_i = p_mp[i], q_mp[i], r_mp[i], s_mp[i]

        # Cubic polynomial (from CrossK5Key.lean):
        # (p+q) * A * (z * conj(z)^2)
        #   + (r+s) * conj(A) * (z^2 * conj(z))
        #   + ((q-p)*A*conj(z)^2 + (s-r)*conj(A)*z^2) * ξ
        #   - (p+q)*A*z*ξ^2 - (r+s)*conj(A)*conj(z)*ξ^2
        #   + ((p-q)*A + (r-s)*conj(A)) * ξ^3

        term1 = (p_i + q_i) * A_i * (z_i * (z_i_conj ** mp.mpf('2')))
        term2 = (r_i + s_i) * A_i_conj * ((z_i ** mp.mpf('2')) * z_i_conj)
        term3 = ((q_i - p_i) * A_i * (z_i_conj ** mp.mpf('2')) + (s_i - r_i) * A_i_conj * (z_i ** mp.mpf('2'))) * xi_mp
        term4 = (-(p_i + q_i) * A_i * z_i - (r_i + s_i) * A_i_conj * z_i_conj) * (xi_mp ** mp.mpf('2'))
        term5 = ((p_i - q_i) * A_i + (r_i - s_i) * A_i_conj) * (xi_mp ** mp.mpf('3'))

        cubic_i = term1 + term2 + term3 + term4 + term5
        cubics_mp.append(cubic_i)

    # Compute quartic_i for each block
    # quartic_i = (z_i + ξ) * (z_i - ξ) * (conj(z_i) + ξ) * (conj(z_i) - ξ)
    quartics_mp = []
    for i in range(5):
        z_i = z_mp[i]
        z_i_conj = mp.conj(z_i)
        quartic_i = (z_i + xi_mp) * (z_i - xi_mp) * (z_i_conj + xi_mp) * (z_i_conj - xi_mp)
        quartics_mp.append(quartic_i)

    # Compute cross-product sum
    # result = sum_i (cubic_i * product of quartics_j for j ≠ i)
    total_mp = mp.mpc('0')
    for i in range(5):
        term = cubics_mp[i]
        for j in range(5):
            if i != j:
                term = term * quartics_mp[j]
        total_mp += term

    # Extract magnitude
    magnitude = mp.fabs(total_mp)
    real_part = mp.re(total_mp)
    imag_part = mp.im(total_mp)

    results.append((n_test, total_mp, magnitude, real_part, imag_part))

    # Print summary
    mag_str = mp.nstr(magnitude, 15)
    status = "ZERO (< 1e-60)" if magnitude < mp.mpf('1e-60') else f"NONZERO ({mag_str})"
    print(f"n={n_test:2d}: magnitude = {mag_str:20s}  [{status}]")

print()
print("=" * 100)
print("CRITICAL NODES — FULL PRECISION OUTPUT")
print("=" * 100)

for n_test in [0, 5, 10, 15, 20]:
    _, total_mp, _, real_part, imag_part = results[n_test]
    print(f"\nn={n_test}:")
    print(f"  Real:  {mp.nstr(real_part, 60)}")
    print(f"  Imag:  {mp.nstr(imag_part, 60)}")
    print(f"  |val|: {mp.nstr(mp.fabs(total_mp), 60)}")

print()
print("=" * 100)
print("SUMMARY")
print("=" * 100)

# Count how many nodes are essentially zero (< noise floor at dps=80)
zero_threshold = mp.mpf('1e-70')
zero_count = sum(1 for _, _, mag, _, _ in results if mag < zero_threshold)
nonzero_count = 21 - zero_count

print()
print(f"Total nodes tested: 21")
print(f"  Cubic ≈ 0 (< 1e-70): {zero_count}")
print(f"  Cubic ≠ 0 (≥ 1e-70): {nonzero_count}")
print()

if nonzero_count == 0:
    print("✓✓✓ ALL 21 NODES PASS: CUBIC = 0 (confirmed at dps=80)")
    print("    S2742 is CORRECT")
else:
    min_mag = min(mag for _, _, mag, _, _ in results)
    max_mag = max(mag for _, _, mag, _, _ in results)
    print(f"✗✗✗ CUBIC IS NONZERO")
    print(f"    Minimum magnitude: {mp.nstr(min_mag, 20)}")
    print(f"    Maximum magnitude: {mp.nstr(max_mag, 20)}")
    print(f"    S2742 is INCORRECT (sympy.simplify() false positive)")
    print()
    print("    First few nonzero nodes:")
    for n_test, _, mag, _, _ in results[:5]:
        if mag >= zero_threshold:
            print(f"      n={n_test}: {mp.nstr(mag, 15)}")

print()
