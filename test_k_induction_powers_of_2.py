"""
Test whether the K-induction cubic/quartic polynomial identity holds on the
"Powers of 2" survivor family: z = [1, 2, 4, 8, 16, 32] + 2.3j

Strategy: Extract the polynomial from CrossK6Key.lean and evaluate at high
precision (mpmath dps=100) on the Powers of 2 z-values.

Result interpretation:
- If sum = 0 at dps=100, the identity holds on this restricted family (partial positive result)
- If sum ≠ 0 even at dps=100, K-induction still fails even on high-Im/independent families
"""

import mpmath as mp
import numpy as np

def evaluate_cubic(n_val, N, L, z_vals_one_block, p, q, r, s, A, A_conj):
    """
    Evaluate the cubic expression for one block.

    cubic_i = (p + q) * A * z_i * conj(z_i)^2
            + (r + s) * conj(A) * z_i^2 * conj(z_i)
            + (q - p) * A * conj(z_i)^2 * (2π(n - N) / L)
            + (s - r) * conj(A) * z_i^2 * (2π(n - N) / L)
            - (p + q) * A * z_i * (2π(n - N) / L)^2
            - (r + s) * conj(A) * conj(z_i) * (2π(n - N) / L)^2
            + (p - q) * A * (2π(n - N) / L)^3
            + (r - s) * conj(A) * (2π(n - N) / L)^3
    """
    z_i = z_vals_one_block
    z_i_conj = mp.conj(z_i)

    xi = mp.mpf(2) * mp.pi * (mp.mpf(n_val) - mp.mpf(N)) / mp.mpf(L)  # 2π(n - N) / L
    xi2 = xi ** 2
    xi3 = xi ** 3

    term1 = (p + q) * A * z_i * (z_i_conj ** 2)
    term2 = (r + s) * A_conj * (z_i ** 2) * z_i_conj
    term3 = (q - p) * A * (z_i_conj ** 2) * xi
    term4 = (s - r) * A_conj * (z_i ** 2) * xi
    term5 = -(p + q) * A * z_i * xi2
    term6 = -(r + s) * A_conj * z_i_conj * xi2
    term7 = (p - q) * A * xi3
    term8 = (r - s) * A_conj * xi3

    return term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8

def evaluate_quartic(n_val, N, L, z_vals_one_block):
    """
    Evaluate the quartic expression for one block.

    quartic_i = (z_i + 2π(n - N) / L) * (z_i - 2π(n - N) / L)
              * (conj(z_i) + 2π(n - N) / L) * (conj(z_i) - 2π(n - N) / L)
    """
    z_i = z_vals_one_block
    z_i_conj = mp.conj(z_i)
    xi = mp.mpf(2) * mp.pi * (mp.mpf(n_val) - mp.mpf(N)) / mp.mpf(L)

    return (z_i + xi) * (z_i - xi) * (z_i_conj + xi) * (z_i_conj - xi)

if __name__ == "__main__":
    mp.dps = 100  # High precision: 100 decimal places

    # Powers of 2 family: z = [1, 2, 4, 8, 16, 32] + 2.3j
    z_values = [mp.mpc(float(2**k), 2.3) for k in range(6)]

    print("=" * 100)
    print("K-INDUCTION POLYNOMIAL TEST: POWERS OF 2 FAMILY")
    print("=" * 100)
    print(f"\nFamily: z = [1, 2, 4, 8, 16, 32] + 2.3j")
    print(f"Precision: mpmath dps=100")
    print(f"z-values: {[str(z) for z in z_values]}")

    # Use arbitrary test parameters (focus is on the identity structure, not the vPlusLC/vMinusLC coefficients)
    # For simplicity, assume coefficients p_i = q_i = r_i = s_i = 1 for all blocks
    # (These would come from the vPlusLC/vMinusLC linear independence constraint in reality)

    # Use arbitrary test values
    N = 10
    L = mp.mpf(20)
    n_test = 15  # Test at one value of n

    # Arbitrary but non-zero coefficients
    p_vals = [mp.mpc(1, 0)] * 6
    q_vals = [mp.mpc(1, 0)] * 6
    r_vals = [mp.mpc(1, 0)] * 6
    s_vals = [mp.mpc(1, 0)] * 6

    print(f"\nTest parameters: N={N}, L={L}, n={n_test}")
    print(f"Coefficients: p_i = q_i = r_i = s_i = 1 (simplified test)")

    # Compute A_i = sqrt(2) * sin(L * z_i / 2) for each block
    sqrt2 = mp.sqrt(mp.mpf(2))
    A_vals = [sqrt2 * mp.sin(L * z / mp.mpf(2)) for z in z_values]
    A_conj_vals = [mp.conj(A) for A in A_vals]

    # Compute cubics and quartics
    cubics = []
    quartics = []

    for i in range(6):
        cubic = evaluate_cubic(n_test, N, L, z_values[i], p_vals[i], q_vals[i], r_vals[i], s_vals[i], A_vals[i], A_conj_vals[i])
        quartic = evaluate_quartic(n_test, N, L, z_values[i])
        cubics.append(cubic)
        quartics.append(quartic)

    # Compute the sum-of-products: cubic_i * (all other quartics)
    total_sum = mp.mpc(0)

    for i in range(6):
        # Product of all quartics except quartic_i
        product_others = mp.mpc(1)
        for j in range(6):
            if i != j:
                product_others *= quartics[j]

        term = cubics[i] * product_others
        total_sum += term

        if i == 0:
            print(f"\nTermwise evaluation (first term shown in detail):")
            print(f"  cubic[0] = {mp.nstr(cubics[0], 15)}")
            print(f"  quartic[0] = {mp.nstr(quartics[0], 15)}")
            print(f"  product of other quartics = {mp.nstr(product_others, 15)}")
            print(f"  term[0] = {mp.nstr(term, 15)}")

    # Final result
    print(f"\nFinal sum (K-induction polynomial):")
    print(f"  {mp.nstr(total_sum, 30)}")

    abs_sum = mp.fabs(total_sum)
    print(f"\nAbsolute value of sum:")
    print(f"  {mp.nstr(abs_sum, 30)}")

    # Check whether it's close to zero
    threshold = mp.mpf(10) ** (-50)  # Very tight threshold for dps=100

    if abs_sum < threshold:
        print(f"\n*** RESULT: IDENTITY HOLDS (sum ≈ 0 at dps=100) ***")
        print(f"This is a GENUINE PARTIAL POSITIVE RESULT on the Powers-of-2 family.")
        print(f"Conclusion: K-induction is NOT architecturally impossible; it works on this restricted parameter space.")
    else:
        print(f"\n*** RESULT: IDENTITY FAILS (sum ≠ 0 at dps=100) ***")
        print(f"Reason: Even with high-Im and independent-vector structure,")
        print(f"  the K-induction polynomial does not vanish on this family.")
        print(f"  This suggests structural obstructions beyond the Im-threshold phenomenon.")

    print("\n" + "=" * 100)
