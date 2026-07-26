#!/usr/bin/env python3
"""ATTACK B3, second pass: Weil form in POSITION space with POSITIVE kernel.

Key fixes from v1:
1. Use |Ω̌(x)| = e^{−|x|/2}/(1−e^{−2|x|}) (absolute value, now positive)
2. Better matrix square root via SVD (handles near-singular and ill-conditioned)
3. Superpotential extraction via finite-difference: W[i] ≈ (A[i,i+1] − A[i-1,i]) / dx
   (extracts the potential from the near-diagonal slope)
4. Test antisymmetry more carefully: form odd/even decomposition
5. Report near-null eigenvector (should relate to ψ₀)
"""
import numpy as np
import sys
sys.path.insert(0, "operators")
from critical_eigvec_id import prime_powers_up_to

def weil_kernel_positive(x):
    """Ω̌(x) = e^{−|x|/2} / (1 − e^{−2|x|}), positive."""
    ax = np.abs(x)
    if ax < 1e-12:
        return 0.25  # limit as x->0
    num = np.exp(-ax / 2.0)
    den = 1.0 - np.exp(-2.0 * ax)
    return num / den if abs(den) > 1e-12 else 0.0

weil_kernel_vec = np.vectorize(weil_kernel_positive)

def build_weil_matrix_v2(xgrid, prime_powers, c_max):
    """Build Weil matrix with POSITIVE kernel + arch boundary."""
    D = len(xgrid)
    T = np.zeros((D, D))

    # Arch kernel: Ω̌(x_i - x_j), positive
    for i in range(D):
        for j in range(D):
            dx_ij = xgrid[i] - xgrid[j]
            T[i, j] += weil_kernel_positive(dx_ij)

    # Boundary rank-2: mild damping at edges (regularization)
    T[0, 0] += 0.5
    T[-1, -1] += 0.5

    # Prime sector: point masses at ±log(q) with weight Λ(p)p^{-1/2}
    for q, p in prime_powers:
        log_q = np.log(q)
        weight = np.log(p) * (q ** -0.5)

        for sign in [1, -1]:
            target_x = sign * log_q
            idx = np.argmin(np.abs(xgrid - target_x))
            if 0 <= idx < D:
                T[idx, idx] += weight

    return T

def sqrtm_via_svd(M):
    """Compute A such that A² ≈ M via SVD. Handles indefinite/singular matrices."""
    U, s, Vt = np.linalg.svd(M)
    # For real symmetric M, U ≈ V, but use SVD to be robust
    s_sqrt = np.sqrt(np.maximum(s, 1e-14))  # clip small negatives
    return U @ np.diag(s_sqrt) @ Vt

def sqrtm_eig(M):
    """Eigenvalue method: M = P D P^T, then sqrt(M) = P sqrt(D) P^T."""
    w, V = np.linalg.eigh(M)
    w = np.maximum(w, 1e-14)
    return V @ np.diag(np.sqrt(w)) @ V.T

def extract_superpotential_v2(A, xgrid):
    """Extract W(x) via finite differences.

    For discrete A ≈ d/dx + W(x) (one-step FD):
    The matrix elements encode: A[i,i] = potential part + kinetic part
                                 A[i,i±1] = off-diag kinetic coupling

    Rough extraction: W[i] ≈ mean(|A[i,i-1]|, |A[i,i+1]|) / dx [kinetic scale]
    Or better: look at anti-Hermitian part of A (but A should be real symmetric)

    Simpler robust method: W[i] = A[i,i] itself (contains potential)
    """
    D = len(A)
    dx = xgrid[1] - xgrid[0] if D > 1 else 1.0

    diag = np.diag(A)

    # Try to remove kinetic part: estimate from off-diag
    # For a true d/dx, off-diag is O(1/dx). Compute mean off-diag scale.
    offdiag_scale = 0.0
    count = 0
    for i in range(D - 1):
        if i < D - 1:
            offdiag_scale += abs(A[i, i + 1])
            count += 1
    if count > 0:
        offdiag_scale /= count
        offdiag_scale /= dx  # normalize by grid spacing

    # Estimate: W = diag - C * (kinetic_scale * dx²)
    # This is crude; a better approach would fit from (A²)[i,i]
    # But for now, report diag directly + note caveat

    return diag, offdiag_scale

def decompose_even_odd(vals, xgrid):
    """Decompose vals(x) into even + odd parts."""
    D = len(vals)
    even = np.zeros(D)
    odd = np.zeros(D)
    for i in range(D):
        mirror_i = D - 1 - i
        if mirror_i >= 0 and mirror_i < D:
            val_mirror = vals[mirror_i]
            even[i] = (vals[i] + val_mirror) / 2.0
            odd[i] = (vals[i] - val_mirror) / 2.0
    return even, odd

def analyze_superpotential_v2(W, xgrid):
    """Analyze W for arithmetic structure."""
    # Decompose into even + odd
    W_even, W_odd = decompose_even_odd(W, xgrid)

    stats = {
        'mean': float(np.mean(W)),
        'std': float(np.std(W)),
        'min': float(np.min(W)),
        'max': float(np.max(W)),
        'even_norm': float(np.linalg.norm(W_even)),
        'odd_norm': float(np.linalg.norm(W_odd)),
        'odd_to_total': float(np.linalg.norm(W_odd) / (np.linalg.norm(W) + 1e-12)),
    }

    if len(W) > 2:
        d2W = np.diff(W, 2)
        stats['rms_curvature'] = float(np.sqrt(np.mean(d2W ** 2)))
        stats['smoothness'] = float(np.sqrt(np.mean(d2W ** 2)) / (np.mean(np.abs(W)) + 1e-12))
    else:
        stats['rms_curvature'] = np.nan
        stats['smoothness'] = np.nan

    return stats

def main():
    c_max = 14
    N_grid = 128  # higher res for better discretization
    x_range = np.log(c_max)

    xgrid = np.linspace(-x_range, x_range, N_grid)
    prime_powers = prime_powers_up_to(c_max)
    dx = xgrid[1] - xgrid[0]

    print("=== Weil Position-Space SUSY v2: POSITIVE kernel ===\n")
    print(f"Config: c_max={c_max}, N={N_grid}, x ∈ [{xgrid[0]:.3f}, {xgrid[-1]:.3f}], dx={dx:.4f}\n")

    # Build matrix
    print("Building Weil matrix (positive kernel + prime masses + boundary)...")
    T = build_weil_matrix_v2(xgrid, prime_powers, c_max)

    sym_err = np.linalg.norm(T - T.T) / (np.linalg.norm(T) + 1e-12)
    print(f"  Symmetry: {sym_err:.2e}")

    eigs_T = np.linalg.eigvalsh(T)
    n_neg = np.sum(eigs_T < -1e-10)
    n_zero = np.sum(np.abs(eigs_T) < 1e-10)
    print(f"  Eigenvalue range: [{eigs_T.min():.4f}, {eigs_T.max():.4f}]")
    print(f"  Negative: {n_neg}, near-zero: {n_zero}")
    print()

    # Report near-null eigenvector (if it exists)
    min_eig_idx = np.argmin(np.abs(eigs_T))
    min_eig_val = eigs_T[min_eig_idx]
    print(f"  Smallest |eigenvalue|: {min_eig_val:.6f} at index {min_eig_idx}")
    if abs(min_eig_val) < 0.1:
        # Eigenvector is near-null, might be ψ₀ structure
        _, Q = np.linalg.eigh(T)
        psi_min = Q[:, min_eig_idx]
        print(f"  Near-null eigenvector norm: {np.linalg.norm(psi_min):.6f}")
        print(f"  Decay pattern (|ψ₀| at edges): {abs(psi_min[0]):.4f}, center: {abs(psi_min[N_grid//2]):.4f}, edge: {abs(psi_min[-1]):.4f}")

    print()

    # Matrix square root
    print("Computing matrix square root A via eigendecomposition...")
    A = sqrtm_eig(T)

    A2 = A @ A
    err_sq = np.linalg.norm(A2 - T) / (np.linalg.norm(T) + 1e-12)
    print(f"  Reconstruction error ‖A² - T‖/‖T‖: {err_sq:.4f}")
    sym_A = np.linalg.norm(A - A.T) / (np.linalg.norm(A) + 1e-12)
    print(f"  A symmetry error: {sym_A:.2e}")
    print()

    # Extract superpotential
    print("Extracting superpotential W(x)...")
    W, offdiag_scale = extract_superpotential_v2(A, xgrid)

    print(f"  Kinetic off-diagonal scale (normalized): {offdiag_scale:.4f}")
    print(f"  W(x) sample values:")
    for idx in [0, N_grid//4, N_grid//2, 3*N_grid//4, N_grid-1]:
        print(f"    W({xgrid[idx]:7.3f}) = {W[idx]:10.6f}")

    print()

    # Analyze
    stats = analyze_superpotential_v2(W, xgrid)
    print(f"W(x) statistics:")
    for k, v in stats.items():
        if isinstance(v, float):
            if np.isnan(v):
                print(f"  {k}: NaN")
            else:
                print(f"  {k}: {v:.6f}")

    print()

    # Decompose W into even + odd
    W_even, W_odd = decompose_even_odd(W, xgrid)
    print(f"Even/odd decomposition:")
    print(f"  ‖W_even‖ / ‖W_total‖ = {np.linalg.norm(W_even)/np.linalg.norm(W):.4f}")
    print(f"  ‖W_odd‖ / ‖W_total‖ = {np.linalg.norm(W_odd)/np.linalg.norm(W):.4f}")
    print(f"  (ζ′/ζ is odd, so odd_norm should be significant if clean structure)")

    print()
    print("=== Arithmetic structure test ===")

    # Check: is W close to a simple transcendental?
    # ζ′/ζ has poles at non-trivial zeros (and at s=1 in the complex plane)
    # In log-space around 0, ζ′(1/2+it)/ ζ(1/2+it) is smooth

    # Numerical check: polynomial fit
    from numpy.polynomial import Polynomial
    try:
        # Fit odd part to polynomial
        x_pos = xgrid[N_grid//2:]
        y_odd = W_odd[N_grid//2:]
        if len(x_pos) > 2 and np.std(y_odd) > 1e-8:
            p_odd = Polynomial.fit(x_pos, y_odd, 3)
            residuals = y_odd - p_odd(x_pos)
            rms_residual = np.sqrt(np.mean(residuals ** 2))
            print(f"  Fit W_odd to degree-3 poly on [0, x_max]:")
            print(f"    RMS residual: {rms_residual:.6f}")
            print(f"    Poly coeff (const, lin, quad, cubic): {p_odd.coef[:4]}")
    except Exception as e:
        print(f"  Poly fit failed: {e}")

    print()
    print("=== Discretization Caveats ===")
    print(f"• Grid spacing: dx = {dx:.4f}")
    print(f"• Weil kernel sampled at {N_grid} points (naive quadrature)")
    print(f"• Prime masses quantized to nearest grid point")
    print(f"• Matrix sqrt via eigendecomposition (stable but masks high-freq artifacts)")
    print(f"• Superpotential extraction: diagonal of A (contains kinetic + potential)")
    print(f"  True W requires deconvolving d/dx part → not done rigorously here")
    print(f"• Antisymmetry/odd-part test is post-hoc, not enforced during construction")

    print("\n" + "="*60)

    # Save
    import json
    results = {
        'config': {'c_max': c_max, 'N': N_grid, 'dx': float(dx)},
        'xgrid': xgrid.tolist(),
        'W': W.tolist(),
        'W_even': W_even.tolist(),
        'W_odd': W_odd.tolist(),
        'A_diag': np.diag(A).tolist(),
        'eigs_T': eigs_T.tolist(),
        'stats': stats,
    }
    outfile = '/tmp/weil_position_space_v2_results.json'
    with open(outfile, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved: {outfile}")

if __name__ == '__main__':
    main()
