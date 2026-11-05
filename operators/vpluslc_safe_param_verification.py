"""
Verify whether Hermite-roots and unit-circle parameters avoid Im-threshold collapse (S2718).

S2718 characterized K=6 collapse boundary: Im(z_i) proximity matters, threshold ~1.92.
S2729 found that Hermite roots and unit circle avoid Vandermonde collapse.
Question: Do these "safe" parameter families also avoid the Im-threshold mechanism?

If yes: Restricted K-induction proof on these families becomes tractable.
"""
import numpy as np

def vPlusLC_simplified(z, L, N, n):
    """
    Simplified vPlusLC evaluation for Im-threshold check.
    Formula: sqrt(2)*sin(L*z/2) / (L*(z + 2π(n-N)/L))

    For characterizing threshold, we focus on: does the denominator get near zero?
    """
    denom = L * (z + 2*np.pi*(n - N) / L)
    if abs(denom) < 1e-10:
        return float('nan')
    return np.sqrt(2) * np.sin(L*z/2) / denom

def check_im_threshold(z_values, K, L_param=None, threshold=1.92):
    """
    Check whether z_values violate the Im-threshold criterion (S2718).

    S2718 hypothesis: collapse occurs when Im(z_i) < threshold ≈ 1.92
    """
    if L_param is None:
        L_param = 2*np.pi*K / (K+1)

    N = K  # test parameter

    # Check Im-threshold for each z_i
    im_parts = np.array([z.imag if isinstance(z, complex) else 0 for z in z_values])
    min_im = np.min(np.abs(im_parts))  # minimum imaginary part (magnitude)

    violates = min_im < threshold

    return min_im, violates

def evaluate_hermite_roots(K, threshold=1.92):
    """Get Hermite orthogonal polynomial roots, check Im-threshold."""
    from scipy.special import roots_hermite

    h_roots, _ = roots_hermite(K)
    # Scale to unit region
    h_roots_scaled = h_roots / max(abs(h_roots)) * 0.5

    # These should be REAL, so Im=0, which WOULD violate threshold
    min_im, violates = check_im_threshold(h_roots_scaled, K, threshold=threshold)

    return h_roots_scaled, min_im, violates

def evaluate_unit_circle(K, threshold=1.92):
    """Get unit circle roots (roots of unity), check Im-threshold."""
    circ_roots = np.array([np.exp(2j*np.pi*k/K) for k in range(K)])

    # On unit circle: |Im(z)| = sin(2πk/K), which varies
    min_im, violates = check_im_threshold(circ_roots, K, threshold=threshold)

    return circ_roots, min_im, violates

def evaluate_standard_params(K, threshold=1.92):
    """Standard params from S2715 / S2729."""
    standard_z = np.array([0.3+0.2j, 0.4+0.15j, 0.5+0.1j, 0.6+0.05j, 0.7+0.0j, 0.8-0.05j])

    min_im, violates = check_im_threshold(standard_z, K, threshold=threshold)

    return standard_z, min_im, violates

if __name__ == "__main__":
    K = 6
    threshold = 1.92  # S2718 characterization

    print(f"Im-Threshold Verification (K={K}, threshold={threshold})")
    print("=" * 100)
    print(f"{'Parameter Family':40s} | {'min|Im(z)|':15s} | {'Violates?':10s} | {'Implication':30s}")
    print("-" * 100)

    # Test 1: Hermite roots (REAL)
    h_roots, h_min_im, h_violates = evaluate_hermite_roots(K, threshold)
    status_h = "VIOLATES" if h_violates else "SAFE"
    print(f"{'Hermite roots (scaled)':40s} | {h_min_im:15.4f} | {status_h:10s} | {'All real → Im=0 → UNSAFE':30s}")

    # Test 2: Unit circle
    circ_roots, circ_min_im, circ_violates = evaluate_unit_circle(K, threshold)
    status_circ = "VIOLATES" if circ_violates else "SAFE"
    print(f"{'Unit circle (exp(2πik/K))':40s} | {circ_min_im:15.4f} | {status_circ:10s} | {'|Im| varies with k':30s}")

    # Test 3: Standard params
    std_params, std_min_im, std_violates = evaluate_standard_params(K, threshold)
    status_std = "VIOLATES" if std_violates else "SAFE"
    print(f"{'Standard (S2715)':40s} | {std_min_im:15.4f} | {status_std:10s} | {'Mixed Im parts':30s}")

    print("=" * 100)

    # Analysis
    print(f"\nAnalysis:")
    print(f"1. Hermite roots: all REAL (Im=0) → violate threshold (0 < {threshold})")
    print(f"   → CANNOT use Hermite roots to avoid Im-threshold collapse")
    print(f"2. Unit circle: |Im(z)| = sin(2πk/K) for k=0..{K-1}")
    circ_ims = [abs(np.sin(2*np.pi*k/K)) for k in range(K)]
    print(f"   → min|Im| = {min(circ_ims):.4f}" + (f" < {threshold}" if min(circ_ims) < threshold else f" >= {threshold}"))
    print(f"   → Unit circle: {'VIOLATES' if min(circ_ims) < threshold else 'SAFE from'} Im-threshold")
    print(f"3. Standard params: min|Im| = {std_min_im:.4f}" + (f" < {threshold}" if std_min_im < threshold else f" >= {threshold}"))

    print(f"\nConclusion:")
    if h_violates and std_violates:
        print(f"  Both Hermite and standard params violate Im-threshold.")
        if circ_violates:
            print(f"  Even unit circle violates (min|Im|={circ_min_im:.4f} < {threshold}).")
            print(f"  → Im-threshold mechanism may be unavoidable; no parameter family escapes it.")
        else:
            print(f"  BUT unit circle does NOT violate (min|Im|={circ_min_im:.4f} >= {threshold}).")
            print(f"  → Unit circle is SAFE; restricted K-induction on unit circle may avoid collapse.")
    else:
        print(f"  At least one family avoids Im-threshold collapse.")
