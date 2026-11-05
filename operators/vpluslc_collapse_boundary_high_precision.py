"""
HIGH-PRECISION determinant check for vPlusLC collapse boundary.

The original vpluslc_collapse_boundary.py used numpy/float64 determinants.
For magnitudes ranging 1e-63 to 1e+2, float64 is INADEQUATE.

This script re-computes key determinants at dps=50, 100 using mpmath
to verify whether the float64 results are real or artifacts.
"""
import mpmath as mp
import numpy as np

def vandermonde_matrix_mpmath(z_values, K, dps=50):
    """Build Vandermonde matrix V[i,j] = z_i^j using mpmath at arbitrary precision."""
    mp.dps = dps
    z_mp = [mp.mpc(complex(z)) for z in z_values]

    V = mp.matrix(K, K)
    for i in range(K):
        for j in range(K):
            V[i, j] = z_mp[i] ** j

    return V

def determinant_comparison(z_values, K, label=""):
    """Compute determinant at float64 and at dps=50, 100 for comparison."""
    print(f"\n{label}")
    print(f"  z_values: {[f'{z:.2f}' for z in z_values]}")

    # Float64 version
    V_float = np.vstack([np.array(z_values)**j for j in range(K)])
    det_float = np.linalg.det(V_float)
    print(f"  Float64 det: {abs(det_float):.4e}")

    # High-precision versions
    for dps in [50, 100]:
        V_mp = vandermonde_matrix_mpmath(z_values, K, dps=dps)
        det_mp = mp.det(V_mp)
        det_mp_abs = float(abs(det_mp))
        print(f"  dps={dps}  det: {det_mp_abs:.4e}")

    return det_float

if __name__ == "__main__":
    K = 6

    print("HIGH-PRECISION DETERMINANT RE-CHECK FOR VPLUSLC COLLAPSE BOUNDARY")
    print("=" * 100)
    print("Testing key parameter families from S2729/S2730 at high precision.")
    print("If float64 and mpmath results differ significantly, float64 is unreliable.")

    # Test 1: Standard parameters (float64 gave 1.84e-10)
    standard_z = np.array([0.3+0.2j, 0.4+0.15j, 0.5+0.1j, 0.6+0.05j, 0.7+0.0j, 0.8-0.05j])
    determinant_comparison(standard_z, K, "TEST 1: Standard params (S2715)")

    # Test 2: Hermite roots (float64 gave 9.16e-07)
    from scipy.special import roots_hermite
    h_roots, _ = roots_hermite(K)
    h_roots_scaled = h_roots / max(abs(h_roots)) * 0.5
    determinant_comparison(h_roots_scaled, K, "TEST 2: Hermite roots (S2729)")

    # Test 3: Unit circle (float64 gave 2.16e+02)
    circ_z = np.array([np.exp(2j*np.pi*k/K) for k in range(K)])
    determinant_comparison(circ_z, K, "TEST 3: Unit circle (S2729)")

    # Test 4: Hermite * 0.5 + 3j (float64 gave 9.16e-07)
    h_shifted = h_roots / max(abs(h_roots)) * 0.5 + 3j
    determinant_comparison(h_shifted, K, "TEST 4: Hermite * 0.5 + 3j (S2730 high-Im)")

    # Test 5: Vertical line Im=5 (float64 gave 7.35e-08)
    vert_z = np.array([j/K + 5j for j in range(1, K+1)])
    determinant_comparison(vert_z, K, "TEST 5: Vertical Im=5 (S2730 high-Im)")

    # Test 6: Dense Im=20 (float64 gave 1.02e-21)
    dense_z = np.array([0.5*np.cos(2*np.pi*k/K) + 1j*20 for k in range(K)])
    determinant_comparison(dense_z, K, "TEST 6: Dense Im=20 (S2730 high-Im)")

    print("\n" + "=" * 100)
    print("INTERPRETATION:")
    print("If float64 and mpmath results agree within 1-2 orders of magnitude,")
    print("  float64 is marginally acceptable (results are real, not pure artifacts).")
    print("If they differ by >2-3 orders of magnitude or show opposite magnitude trends,")
    print("  float64 is UNRELIABLE and S2729/S2730 conclusions are SUSPECT.")
