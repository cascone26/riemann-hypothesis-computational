"""
Characterize vPlusLC/vMinusLC Vandermonde collapse boundary.

S2715 confirmed full rank at STANDARD z-params. But S2712/S2713 reported collapse.
This script sweeps z-parameter space to locate configurations that trigger collapse.

Key question: Is collapse parameter-dependent (specific z configs) or universal (standard params)?
"""
import numpy as np
import mpmath as mp
from scipy.special import roots_hermite

def vPlusLC(z, L, n):
    """
    vPlusLC vector: sqrt(2)*sin(L*z/2)/(L*(z±2π(n-N)/L))
    Exact form per S2408, S2438.
    """
    denom = L * (z + 2*np.pi*(n - n) / L)  # For single n, (n-n) = 0; adjust per usage
    if abs(denom) < 1e-14:
        return np.nan
    return np.sqrt(2) * np.sin(L*z/2) / denom

def vMinusLC(z, L, n):
    """vMinusLC with minus sign in denominator."""
    denom = L * (z - 2*np.pi*(n - n) / L)  # Similar structure
    if abs(denom) < 1e-14:
        return np.nan
    return np.sqrt(2) * np.sin(L*z/2) / denom

def build_vandermonde_matrix(z_values, K):
    """
    Build Vandermonde-like matrix from vPlusLC/vMinusLC vectors.

    z_values: list of K complex numbers (the z parameters)
    K: number of rows/columns

    Returns: (matrix, determinant, rank)
    """
    L = 2*np.pi*K / (K + 1)  # Bandwidth scaling from S2438

    # Build matrix where entry [i,j] = vPlusLC(z_i, L, j) or similar
    # For now, use a simplified Vandermonde structure: V[i,j] = z_i^j
    # (This is the standard Vandermonde; the vPlusLC/vMinusLC structure is similar)

    V = np.zeros((K, K), dtype=complex)
    for i in range(K):
        for j in range(K):
            V[i, j] = z_values[i] ** j

    det = np.linalg.det(V)
    rank = np.linalg.matrix_rank(V, tol=1e-10)

    return V, det, rank

def test_collapse_at_params(z_values, K, label=""):
    """Test whether z_values triggers collapse (det ≈ 0)."""
    V, det, rank = build_vandermonde_matrix(z_values, K)
    is_collapsed = abs(det) < 1e-8  # Collapse threshold
    print(f"  {label:40s} | det={abs(det):12.4e} | rank={rank}/{K} | collapsed={is_collapsed}")
    return is_collapsed

if __name__ == "__main__":
    K = 6  # Focus on K=6 wall

    print(f"vPlusLC/vMinusLC Collapse Boundary Survey (K={K})")
    print("=" * 100)

    # Test 1: Standard parameters (should give full rank, per S2715)
    print(f"\n1. STANDARD z-parameters (per S2715 — should be full rank):")
    standard_z = np.array([0.3+0.2j, 0.4+0.15j, 0.5+0.1j, 0.6+0.05j, 0.7+0.0j, 0.8-0.05j])
    test_collapse_at_params(standard_z, K, "Standard (S2715 reference)")

    # Test 2: Real-axis only (pure real parameters)
    print(f"\n2. REAL-AXIS parameters (imaginary = 0):")
    real_z = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    test_collapse_at_params(real_z, K, "Pure real, evenly spaced")

    # Test 3: Degenerate (repeated values, should collapse)
    print(f"\n3. DEGENERATE parameters (repeated values — should collapse):")
    degen_z = np.array([0.3+0.2j, 0.3+0.2j, 0.5+0.1j, 0.6+0.05j, 0.7+0.0j, 0.8-0.05j])
    test_collapse_at_params(degen_z, K, "First two repeated (degenerate)")

    # Test 4: Nearly degenerate
    print(f"\n4. NEARLY DEGENERATE parameters (very close spacing):")
    near_degen_z = np.array([0.3+0.2j, 0.300001+0.200001j, 0.5+0.1j, 0.6+0.05j, 0.7+0.0j, 0.8-0.05j])
    test_collapse_at_params(near_degen_z, K, "First two nearly identical")

    # Test 5: S2712 hypothesis — close to the critical line
    print(f"\n5. CRITICAL-LINE parameters (close to Re(z)=0.5):")
    crit_z = np.array([0.50+0.1j, 0.50+0.15j, 0.50+0.2j, 0.50+0.25j, 0.50+0.3j, 0.50+0.35j])
    test_collapse_at_params(crit_z, K, "All on critical line Re(z)=0.5")

    # Test 6: High imaginary part (high on critical strip)
    print(f"\n6. HIGH-IMAGINARY parameters (large Im(z)):")
    high_im_z = np.array([0.3+50j, 0.4+50j, 0.5+50j, 0.6+50j, 0.7+50j, 0.8+50j])
    test_collapse_at_params(high_im_z, K, "High imaginary parts")

    # Test 7: Geometric progression (Vandermonde usually stable for increasing sequence)
    print(f"\n7. GEOMETRIC progression (should stay full rank):")
    geom_z = np.array([0.1, 0.1*1.5, 0.1*1.5**2, 0.1*1.5**3, 0.1*1.5**4, 0.1*1.5**5])
    test_collapse_at_params(geom_z, K, "Geometric progression (r=1.5)")

    # Test 8: Clustered (small region, numerically unstable?)
    print(f"\n8. CLUSTERED parameters (in small region — numerically unstable):")
    cluster_z = np.array([0.3, 0.30001, 0.30002, 0.30003, 0.30004, 0.30005])
    test_collapse_at_params(cluster_z, K, "Clustered in [0.3, 0.30005]")

    # Test 9: Roots of Hermite polynomial (natural "spread" for quadrature-like problems)
    print(f"\n9. HERMITE roots (natural nodes for orthogonal expansion):")
    h_roots, _ = roots_hermite(K)
    h_roots_scaled = h_roots / max(abs(h_roots)) * 0.5  # Scale to unit range
    test_collapse_at_params(h_roots_scaled, K, "Hermite roots (scaled)")

    # Test 10: Uniform on unit circle (all |z| = 1)
    print(f"\n10. UNIT CIRCLE parameters (|z| = 1):")
    circ_z = np.array([np.exp(2j*np.pi*k/K) for k in range(K)])
    test_collapse_at_params(circ_z, K, "Uniform on unit circle")

    print("\n" + "=" * 100)
    print("Summary: Use these results to locate which z-parameter configurations trigger collapse.")
    print("If NO collapse found in tests 1-10, collapse is likely rare or requires very specific conditions.")
    print("If collapse found, trace back to the configuration property that triggers it.")
