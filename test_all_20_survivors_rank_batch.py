"""
Batch test: vPlusLC/vMinusLC rank on ALL 20 survivor families.

After fork's structural error correction:
- vPlusLC/vMinusLC rank ≠ hsum constraint rank (different layers)
- vPlusLC/vMinusLC rank = whether z-values produce independent basis vectors
- hsum constraint rank = solution existence for coefficient system
- K=5 and K=6 both have NON-TRIVIAL hsum solutions (9 and 13 free params respectively)

This batch tests whether ALL 20 survivors have the same vPlusLC/vMinusLC rank pattern (6 = full).
"""

import mpmath as mp
import numpy as np
from scipy.linalg import svd

def vPlusLC(n, L, z):
    if isinstance(n, int):
        n = mp.mpf(n)
    if isinstance(L, float):
        L = mp.mpf(L)
    if isinstance(z, complex):
        z = mp.mpc(z)
    sqrt2 = mp.sqrt(mp.mpf(2))
    arg = L * z / mp.mpf(2)
    sin_term = mp.sin(arg)
    leg_n = mp.legendre(n, z)
    return sqrt2 * sin_term * leg_n

def vMinusLC(n, L, z):
    if isinstance(n, int):
        n = mp.mpf(n)
    if isinstance(L, float):
        L = mp.mpf(L)
    if isinstance(z, complex):
        z = mp.mpc(z)
    sqrt2 = mp.sqrt(mp.mpf(2))
    arg = L * z / mp.mpf(2)
    sin_term = mp.sin(arg)
    cheb_n = mp.chebyu(n, z)
    return sqrt2 * sin_term * cheb_n

def test_vpluslc_rank(family_name, z_values, L, N=10):
    """Test vPlusLC/vMinusLC rank on one family."""
    basis_size = 2 * (N + 1)
    num_z = len(z_values)

    matrix = []
    for n in range(N + 1):
        row_plus = [vPlusLC(n, L, z) for z in z_values]
        matrix.append(row_plus)
        row_minus = [vMinusLC(n, L, z) for z in z_values]
        matrix.append(row_minus)

    matrix_float = np.array([[complex(entry) for entry in row]
                             for row in matrix], dtype=np.complex128)

    U, s, Vh = svd(matrix_float, full_matrices=False)
    rank = np.sum(s > 1e-10)

    return rank, num_z

if __name__ == "__main__":
    mp.dps = 50
    L = mp.mpf(20)

    print("=" * 120)
    print("BATCH RANK TEST: ALL 20 SURVIVOR FAMILIES")
    print("=" * 120)

    # All 20 survivors from adversarial search
    families = [
        ("Asymmetric tight-loose (Im=3)",
         [0.0+3j, 0.01+3j, 0.5+3j, 0.51+3j, 1.5+3j, 1.51+3j]),

        ("Asymmetric wide-cluster (Im=2.5)",
         [0.0+2.5j, 0.1+2.5j, 0.2+2.5j, 1.0+2.5j, 1.1+2.5j, 5.0+2.5j]),

        ("Exponential real spacing (Im=2.5)",
         [np.exp(k/2)+2.5j for k in range(6)]),

        ("Powers of 2 (Im=2.3)",
         [2**k + 2.3j for k in range(6)]),

        ("π multiples (Im=3)",
         [k*np.pi + 3j for k in range(1, 7)]),

        ("φ (golden ratio) spacing (Im=2.5)",
         [(1.618034**k) / 10 + 2.5j for k in range(6)]),

        ("√2 spacing (Im=3.1)",
         [k*np.sqrt(2) + 3.1j for k in range(6)]),

        ("e spacing (Im=2.8)",
         [k*np.e + 2.8j for k in range(6)]),

        ("Two clusters (Im=3.5)",
         [0.0+3.5j, 0.01+3.5j, 0.02+3.5j, 10.0+3.5j, 10.01+3.5j, 10.02+3.5j]),

        ("Symmetric Re, varying Im (Im~2.5)",
         [-2.0+2.4j, -1.0+2.5j, 0.0+2.6j, 1.0+2.5j, 2.0+2.4j, 3.0+2.5j]),

        ("Logarithmic spiral (Im > 2.5)",
         [np.exp(k/2)*np.cos(k) + 1j*(2.5 + k*0.3) for k in range(6)]),

        ("Archimedean spiral (Im > 2)",
         [k*np.cos(k)*0.5 + 1j*(2.0 + k*0.4) for k in range(6)]),

        ("Pseudo-random lcg (Im=3)",
         [(((k*1103515245 + 12345) % (2**31)) / 2**31) * 10 + 3j for k in range(6)]),

        ("Hash-based spread (Im=2.5)",
         [((hash(f"z{k}") % 1000) / 100) + 2.5j for k in range(6)]),

        ("Mixed Im (threshold boundary)",
         [0.0+1.95j, 0.1+2.0j, 0.2+2.5j, 0.3+3.0j, 0.4+5.0j, 0.5+10.0j]),

        ("Re modulated by Im",
         [np.sin(k)*np.abs(2.5 + k*0.5) + 1j*(2.5 + k*0.5) for k in range(6)]),

        ("Re=Im tangent (Im > 2)",
         [np.tan(k/10) + 1j*(2.0 + k*0.5) for k in range(6)]),

        ("Extreme high Im=10",
         [k/6 + 10j for k in range(6)]),

        ("Extreme high Im=50",
         [k/6 + 50j for k in range(6)]),

        ("Extreme high Im=100",
         [k/6 + 100j for k in range(6)]),
    ]

    print(f"\nTesting {len(families)} families")
    print(f"{'Family Name':50s} | {'Rank':5s} | {'Result':30s}")
    print("-" * 120)

    full_rank_count = 0
    rank_deficient_count = 0

    for name, z_vals in families:
        # Convert to mpmath
        z_mp = [mp.mpc(complex(z)) for z in z_vals]

        rank, num_z = test_vpluslc_rank(name, z_mp, L)

        if rank == num_z:
            result = "FULL RANK (independent vectors)"
            full_rank_count += 1
        else:
            result = f"RANK-DEFICIENT ({rank}/{num_z})"
            rank_deficient_count += 1

        print(f"{name:50s} | {rank:5d} | {result:30s}")

    print("-" * 120)
    print(f"\nSUMMARY:")
    print(f"  Full rank (linear independence): {full_rank_count} / {len(families)}")
    print(f"  Rank-deficient (linear dependence): {rank_deficient_count} / {len(families)}")

    if full_rank_count == len(families):
        print(f"\n*** UNIVERSAL PATTERN: ALL 20 survivors have full-rank vPlusLC/vMinusLC vectors ***")
        print(f"This means all 20 families produce LINEARLY INDEPENDENT basis vectors.")
        print(f"However: This does NOT prevent non-trivial hsum solutions (as shown by K=5/K=6 comparison).")
    elif rank_deficient_count > 0:
        print(f"\n*** MIXED PATTERN: {rank_deficient_count} families are rank-deficient ***")
        print(f"This means some survivors have linearly dependent vPlusLC/vMinusLC vectors.")
        print(f"Next: Test hsum constraint rank on the rank-deficient families to see if any admit K-induction.")

    print("=" * 120)
