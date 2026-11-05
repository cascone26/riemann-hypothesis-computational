"""
Search for parameter families that STRICTLY avoid Im-threshold (all Im(z) > threshold).

Previous test showed: Hermite roots, unit circle, standard params all have min|Im(z)|=0.
This script tests MODIFIED versions that avoid the real axis entirely.
"""
import numpy as np

def check_im_threshold_strict(z_values, K, threshold=1.92):
    """Check whether ALL components satisfy |Im(z_i)| > threshold."""
    im_parts = np.array([z.imag if isinstance(z, complex) else 0 for z in z_values])
    min_im = np.min(np.abs(im_parts))
    violates = min_im < threshold
    max_im = np.max(np.abs(im_parts))
    return min_im, max_im, violates

def test_family(name, z_values, K, threshold=1.92):
    """Test one parameter family."""
    det_v = np.linalg.det(np.vstack([z_values**j for j in range(K)]))
    min_im, max_im, violates = check_im_threshold_strict(z_values, K, threshold)
    status = "VIOLATES" if violates else "SAFE"
    print(f"{name:45s} | Im ∈ [{min_im:6.4f}, {max_im:6.4f}] | det={abs(det_v):12.4e} | {status:10s}")
    return violates, abs(det_v)

if __name__ == "__main__":
    K = 6
    threshold = 1.92

    print(f"Strictly Safe Parameter Family Search (K={K}, threshold={threshold})")
    print("=" * 130)
    print(f"{'Parameter Family':45s} | {'Im Range':20s} | {'det(V)':15s} | {'Status':10s}")
    print("-" * 130)

    # Test 1: Unit circle shifted upward
    circ_shifted = np.array([np.exp(2j*np.pi*k/K) + 2j for k in range(K)])
    test_family("Unit circle + 2j (shifted up)", circ_shifted, K, threshold)

    # Test 2: Unit circle scaled and shifted
    circ_scaled_shifted = np.array([np.exp(2j*np.pi*k/K) * 0.5 + 3j for k in range(K)])
    test_family("Unit circle * 0.5 + 3j (scaled & shifted)", circ_scaled_shifted, K, threshold)

    # Test 3: Hermite roots shifted upward
    from scipy.special import roots_hermite
    h_roots, _ = roots_hermite(K)
    h_roots_shifted = h_roots / max(abs(h_roots)) * 0.5 + 3j  # pure imaginary shift
    test_family("Hermite roots * 0.5 + 3j (shifted)", h_roots_shifted, K, threshold)

    # Test 4: Vertical line in upper half-plane (Im > 3)
    vert_line = np.array([j/K + 3.5j for j in range(1, K+1)])
    test_family("Vertical line: Re ∈ [1/6, 1], Im=3.5", vert_line, K, threshold)

    # Test 5: Upper half-plane circle (Im > 0 always)
    uhp_circle = np.array([np.cos(2*np.pi*k/K) + 1j*(2+np.sin(2*np.pi*k/K)) for k in range(K)])
    test_family("UHP circle: center at 2i, radius 1", uhp_circle, K, threshold)

    # Test 6: High imaginary, low real spread
    high_im = np.array([j/K + 5j for j in range(K)])
    test_family("High Im: Re ∈ [0, 1], Im=5", high_im, K, threshold)

    # Test 7: Ellipse in upper half-plane
    ellipse_uhp = np.array([0.6*np.cos(2*np.pi*k/K) + 1j*(2+0.8*np.sin(2*np.pi*k/K)) for k in range(K)])
    test_family("Ellipse in UHP: (a,b)=(0.6,0.8), center=2i", ellipse_uhp, K, threshold)

    # Test 8: Extreme high Im (Im > 10)
    extreme_im = np.array([0.1*k + 10j for k in range(K)])
    test_family("Extreme high Im: Im=10", extreme_im, K, threshold)

    # Test 9: Circularly arranged but in upper half-plane
    arc_uhp = np.array([np.exp(1j*(np.pi/2 + np.pi*k/K)) for k in range(K)])  # arc from i to -i via -1
    test_family("Arc in UHP: Im ∈ [0.5, 1]", arc_uhp, K, threshold)

    # Test 10: Dense packing in high-Im region
    dense_high = np.array([0.5*np.cos(2*np.pi*k/K) + 1j*20 for k in range(K)])
    test_family("Dense at high Im: Im=20", dense_high, K, threshold)

    print("=" * 130)
    print(f"\nConclusion:")
    print(f"Search for parameter families with ALL Im(z_i) > {threshold}.")
    print(f"If found: restricted K-induction proof on such families may work.")
    print(f"If none found: Im-threshold mechanism is FUNDAMENTAL, not parameter-escapable.")
