#!/usr/bin/env python3
"""Match observed coefficients to Jacobi theta / Eisenstein q-series."""
import mpmath as mp
mp.mp.dps=20
# Observed coefficients (magnitude, j=0..12)
a_obs = [mp.mpf(x) for x in [0.4409, 0.7882, 0.4103, 0.1246, 0.02208, 0.003189, 0.0003445,
         3.24e-5, 2.871e-6, 2.215e-7, 1.731e-8, 1.302e-9, 7.792e-11]]
q_guess = mp.mpf(0.154)
print("THETA FUNCTION MATCHER")
print("=" * 70)
print("Target: match observed coefficients to known theta/Eisenstein series")
print(f"q-value from decay analysis: {mp.nstr(q_guess, 6)}")
print("=" * 70)
# Test 1: theta_3(0, q) = 1 + 2*sum_{n=1}^inf q^(n^2)
print("\nTEST 1: theta_3(0, q) = 1 + 2·Σ q^(n²)")
c_theta3 = [mp.mpf(0)]  # placeholder for n=0
for n in range(1, 13):
    c_theta3.append(2 * q_guess**(n**2))
print("n\tθ_3 coeff\t\tObserved a_n\t\tRatio")
print("-" * 70)
for n in range(len(a_obs)):
    if n < len(c_theta3):
        ratio = a_obs[n] / c_theta3[n] if c_theta3[n] > 0 else mp.mpf('nan')
        print(f"{n}\t{mp.nstr(c_theta3[n], 8)}\t\t{mp.nstr(a_obs[n], 8)}\t\t{mp.nstr(ratio, 5)}")
# Test 2: theta_4(0, q) = 1 + 2·sum_{n=1}^inf (-q)^(n^2) = 1 + 2·sum (-1)^(n^2)·q^(n^2)
print("\nTEST 2: theta_4(0, q) = 1 + 2·Σ (-1)^(n²)·q^(n²)")
c_theta4 = [mp.mpf(0)]
for n in range(1, 13):
    c_theta4.append(2 * ((-1)**(n**2)) * q_guess**(n**2))
print("n\tθ_4 coeff\t\tObserved a_n\t\tRatio")
print("-" * 70)
for n in range(len(a_obs)):
    if n < len(c_theta4):
        ratio = a_obs[n] / c_theta4[n] if c_theta4[n] != 0 else mp.mpf('nan')
        print(f"{n}\t{mp.nstr(c_theta4[n], 8)}\t\t{mp.nstr(a_obs[n], 8)}\t\t{mp.nstr(ratio, 5)}")
# Test 3: E_4(q) = 1 + 240·sum_{n=1}^inf n^3·q^n / (1 - q^n) [Eisenstein weight 4]
# Simpler: E_4(q) ~ sum a_n q^n with known coefficients
print("\nTEST 3: Eisenstein E_4(q) [normalized]")
# E_4 has the form 1 + 240·σ_3(n)·q^n where σ_3(n) = sum of cubes of divisors
e4_coeffs = []
for n in range(13):
    if n == 0:
        e4_coeffs.append(mp.mpf(1))
    else:
        sigma3 = sum(d**3 for d in range(1, n+1) if n % d == 0)
        e4_coeffs.append(mp.mpf(240) * sigma3 * q_guess**n)
print("n\tE_4(q) coeff\t\tObserved a_n\t\tRatio")
print("-" * 70)
for n in range(len(a_obs)):
    if n < len(e4_coeffs):
        ratio = a_obs[n] / e4_coeffs[n] if e4_coeffs[n] > 0 else mp.mpf('nan')
        print(f"{n}\t{mp.nstr(e4_coeffs[n], 8)}\t\t{mp.nstr(a_obs[n], 8)}\t\t{mp.nstr(ratio, 5)}")
print("\n" + "=" * 70)
print("BEST-FIT ANALYSIS:")
print("If all ratios are CONSTANT across n, we have a match (constant = scaling factor)")
print("Theta_3: (check if ratio ~const); Theta_4: (check); E_4: (check)")
print("=" * 70)
