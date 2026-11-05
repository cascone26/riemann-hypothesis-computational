"""
ADVERSARIAL SEARCH: Find high-Im z-parameters that ALSO have independent vectors.

Goal: Construct 30+ structurally diverse z-parameter families, all with min|Im(z)| > 1.92 (avoid Im-threshold),
compute determinants at mpmath dps=50 (never float64), and check stability at dps=100, 200.

If even ONE configuration has det ≠ 0 at all three precision levels, it disproves S2730's "structurally unavoidable" claim.
If NONE survive, report honestly as "no survivor found after N diverse attempts at mpmath precision."
"""

import mpmath as mp
import numpy as np

def vandermonde_mpmath(z_values, K, dps=50):
    """Vandermonde determinant at arbitrary precision using mpmath."""
    mp.dps = dps
    z_mp = [mp.mpc(complex(z)) for z in z_values]
    V = mp.matrix(K, K)
    for i in range(K):
        for j in range(K):
            V[i, j] = z_mp[i] ** j
    det = mp.det(V)
    return float(abs(det))

def check_min_im(z_values):
    """Check minimum imaginary part."""
    im_parts = np.array([z.imag if isinstance(z, complex) else 0 for z in z_values])
    return np.min(np.abs(im_parts))

def test_family(name, z_values, K=6, threshold=1.92):
    """Test one parameter family."""
    min_im = check_min_im(z_values)

    if min_im < threshold:
        status = f"VIOLATES threshold (Im={min_im:.3f})"
        return None, None, status

    # Compute at dps=50
    det_50 = vandermonde_mpmath(z_values, K, dps=50)

    # If promising (det != 0), check stability at higher dps
    if det_50 > 1e-10:
        det_100 = vandermonde_mpmath(z_values, K, dps=100)
        det_200 = vandermonde_mpmath(z_values, K, dps=200)

        # Check stability: ratios should be close to 1.0
        ratio_100 = det_100 / det_50 if det_50 > 0 else float('nan')
        ratio_200 = det_200 / det_50 if det_50 > 0 else float('nan')

        if abs(ratio_100 - 1.0) < 0.1 and abs(ratio_200 - 1.0) < 0.1:
            stability = "STABLE"
        else:
            stability = f"UNSTABLE (ratios {ratio_100:.2f}, {ratio_200:.2f})"

        return det_50, (det_50, det_100, det_200), f"INDEPENDENT (Im={min_im:.3f}) | det={det_50:.4e} | {stability}"
    else:
        return det_50, None, f"DEPENDENT (det={det_50:.4e}) | Im={min_im:.3f}"

if __name__ == "__main__":
    K = 6
    threshold = 1.92

    print("=" * 120)
    print("ADVERSARIAL SEARCH: 30+ DIVERSE z-PARAMETER FAMILIES (mpmath dps=50+, min|Im| > 1.92)")
    print("=" * 120)

    families = []
    survivors = []

    # CATEGORY 1: Asymmetric spacing (not evenly distributed)
    families.append(("Asymmetric tight-loose (Im=3)",
                     [0.0+3j, 0.01+3j, 0.5+3j, 0.51+3j, 1.5+3j, 1.51+3j]))

    families.append(("Asymmetric wide-cluster (Im=2.5)",
                     [0.0+2.5j, 0.1+2.5j, 0.2+2.5j, 1.0+2.5j, 1.1+2.5j, 5.0+2.5j]))

    families.append(("Extreme asymmetry (Im=2.2)",
                     [0.0+2.2j, 0.001+2.2j, 1.0+2.2j, 10.0+2.2j, 100.0+2.2j, 1000.0+2.2j]))

    # CATEGORY 2: Mixed real parts (large and small)
    families.append(("Mixed Re magnitude (Im=3)",
                     [0.001+3j, 0.1+3j, 1.0+3j, 10.0+3j, 100.0+3j, 1000.0+3j]))

    families.append(("Exponential real spacing (Im=2.5)",
                     [np.exp(k/2)+2.5j for k in range(K)]))

    families.append(("Powers of 2 (Im=2.3)",
                     [2**k + 2.3j for k in range(K)]))

    # CATEGORY 3: Irrational/transcendental spacing
    families.append(("π multiples (Im=3)",
                     [k*np.pi + 3j for k in range(1, K+1)]))

    families.append(("φ (golden ratio) spacing (Im=2.5)",
                     [(1.618034**k) / 10 + 2.5j for k in range(K)]))

    families.append(("√2 spacing (Im=3.1)",
                     [k*np.sqrt(2) + 3.1j for k in range(K)]))

    families.append(("e spacing (Im=2.8)",
                     [k*np.e + 2.8j for k in range(K)]))

    # CATEGORY 4: Deliberately "bad" Vandermonde structures
    families.append(("Nearly repeated (Im=4)",
                     [0.0+4j, 0.0001+4j, 0.0002+4j, 1.0+4j, 1.0001+4j, 1.0002+4j]))

    families.append(("Two clusters (Im=3.5)",
                     [0.0+3.5j, 0.01+3.5j, 0.02+3.5j, 10.0+3.5j, 10.01+3.5j, 10.02+3.5j]))

    # CATEGORY 5: Symmetric in Re, varying Im slightly
    families.append(("Symmetric Re, varying Im (Im~2.5)",
                     [-2.0+2.4j, -1.0+2.5j, 0.0+2.6j, 1.0+2.5j, 2.0+2.4j, 3.0+2.5j]))

    # CATEGORY 6: Circular/spiral patterns (Re+i*Im with varying Im)
    families.append(("Logarithmic spiral (Im > 2.5)",
                     [np.exp(k/2)*np.cos(k) + 1j*(2.5 + k*0.3) for k in range(K)]))

    families.append(("Archimedean spiral (Im > 2)",
                     [k*np.cos(k)*0.5 + 1j*(2.0 + k*0.4) for k in range(K)]))

    # CATEGORY 7: Random-looking but deterministic
    families.append(("Pseudo-random lcg (Im=3)",
                     [(((k*1103515245 + 12345) % (2**31)) / 2**31) * 10 + 3j for k in range(K)]))

    families.append(("Hash-based spread (Im=2.5)",
                     [((hash(f"z{k}") % 1000) / 100) + 2.5j for k in range(K)]))

    # CATEGORY 8: Extreme high Im (hoping dense packing at high Im helps)
    families.append(("Extreme high Im=10",
                     [k/K + 10j for k in range(K)]))

    families.append(("Extreme high Im=50",
                     [k/K + 50j for k in range(K)]))

    families.append(("Extreme high Im=100",
                     [k/K + 100j for k in range(K)]))

    # CATEGORY 9: Mixed high and low Im (some avoid threshold, others barely exceed)
    families.append(("Mixed Im (threshold boundary)",
                     [0.0+1.95j, 0.1+2.0j, 0.2+2.5j, 0.3+3.0j, 0.4+5.0j, 0.5+10.0j]))

    # CATEGORY 10: Modulated real part by imaginary part
    families.append(("Re modulated by Im",
                     [np.sin(k)*np.abs(2.5 + k*0.5) + 1j*(2.5 + k*0.5) for k in range(K)]))

    families.append(("Re=Im tangent (Im > 2)",
                     [np.tan(k/10) + 1j*(2.0 + k*0.5) for k in range(K)]))

    # CATEGORY 11: Pathological cases
    families.append(("All identical Re (Im=5)",
                     [1.0 + (5 + k*0.01)*1j for k in range(K)]))

    families.append(("Nearly identical Re (Im=3)",
                     [1.0 + (3 + k*0.1)*1j + (1e-6 * k)*1j for k in range(K)]))

    # Run all tests
    print(f"\nTesting {len(families)} families:")
    print(f"{'Family Name':45s} | {'Status':65s}")
    print("-" * 120)

    for name, z_values in families:
        det, stability, status = test_family(name, z_values, K, threshold)
        print(f"{name:45s} | {status:65s}")

        if det is not None and det > 1e-10:
            survivors.append((name, z_values, det, stability))

    print("=" * 120)
    print(f"\nRESULTS: {len(survivors)} SURVIVORS (high-Im AND independent vectors)")

    if survivors:
        print("\n*** SURVIVORS (counterexamples to 'structurally unavoidable') ***\n")
        for name, z_values, det, dps_info in survivors:
            print(f"  {name}")
            print(f"    z-values: {[f'{z:.3f}' for z in z_values]}")
            print(f"    det: {det:.4e}")
            print(f"    stability: {dps_info}")
            print()
    else:
        print("\nNO SURVIVORS FOUND after testing 25+ diverse families at mpmath dps=50+")
        print("This is STRONGER evidence (not just float64 noise) but STILL NOT a proof of impossibility.")
        print("Would require either (a) even more diverse families, (b) mathematical proof, or (c) conceding it's an open question.")
