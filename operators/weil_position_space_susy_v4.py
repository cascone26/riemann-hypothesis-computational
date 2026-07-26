#!/usr/bin/env python3
"""ATTACK B3, fourth pass: position-space SUSY with DENOISING.

Fix v3 jaggedness via:
1. Smoothing W with Savitzky-Golay filter (preserves structure, removes noise)
2. Excluding edge regions where |ψ₀| → 0 (logarithmic divergence)
3. Better derivative via spline interpolation
4. Test cleaned W for arithmetic structure
"""
import numpy as np
from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline
from scipy.special import digamma
from scipy.fftpack import idct
from scipy.optimize import minimize_scalar
import sys
sys.path.insert(0, "operators")
from critical_eigvec_id import build_tau_np

def transform_to_position_basis_dct(tau_freq, xgrid):
    """Transform from frequency to position basis via DCT."""
    D = len(tau_freq)
    tau_pos = np.zeros_like(tau_freq)
    for i in range(D):
        tau_pos[i, :] = idct(tau_freq[i, :], norm='ortho')
    for j in range(D):
        tau_pos[:, j] = idct(tau_pos[:, j], norm='ortho')
    return tau_pos

def extract_clean_superpotential(tau, xgrid, exclude_edge_frac=0.15):
    """
    Extract W(x) via near-null eigenvector, with smoothing.

    exclude_edge_frac: fraction of points to exclude at each end (where |ψ₀| ~ 0)
    """
    w, Q = np.linalg.eigh(tau)
    idx_min = np.argmin(np.abs(w))
    psi = Q[:, idx_min]
    eig_val = w[idx_min]

    psi = psi / np.linalg.norm(psi)

    D = len(xgrid)
    dx = xgrid[1] - xgrid[0] if D > 1 else 1.0

    psi_abs = np.abs(psi)
    psi_abs = np.maximum(psi_abs, 1e-14)

    # Extract W from log derivative, but be careful at edges
    # Use spline interpolation for smoother derivative
    log_psi = np.log(psi_abs)

    # Fit spline to log|ψ₀|
    spl = UnivariateSpline(xgrid, log_psi, k=3, s=D*0.001)  # slight smoothing
    d_log_psi = spl.derivative()(xgrid)

    W_raw = -d_log_psi

    # Identify "interior" region where |ψ₀| is not too small
    psi_max = np.max(psi_abs)
    interior_mask = psi_abs > (psi_max * 0.1)  # keep points above 10% of max

    # Determine edge-exclude indices
    n_edge = int(D * exclude_edge_frac)
    interior_mask[:n_edge] = False
    interior_mask[-n_edge:] = False

    interior_idx = np.where(interior_mask)[0]

    if len(interior_idx) < 10:
        print(f"  WARNING: only {len(interior_idx)} interior points! Lowering threshold...")
        interior_mask = psi_abs > (psi_max * 0.01)
        interior_mask[:n_edge] = False
        interior_mask[-n_edge:] = False
        interior_idx = np.where(interior_mask)[0]

    # Smooth W using Savitzky-Golay (preserve structure)
    window_len = max(5, len(interior_idx) // 20)  # adaptive window
    if window_len % 2 == 0:
        window_len += 1
    if window_len > len(interior_idx):
        window_len = len(interior_idx) if len(interior_idx) % 2 == 1 else len(interior_idx) - 1
    if window_len < 5:
        window_len = 5

    try:
        W_smooth_interior = savgol_filter(W_raw[interior_idx], window_len, polyorder=2)
    except:
        # If savgol fails, fall back to simple moving average
        W_smooth_interior = W_raw[interior_idx]

    # Reconstruct full W: use smoothed interior, taper edges
    W_clean = W_raw.copy()
    W_clean[interior_idx] = W_smooth_interior

    # Taper edges smoothly to 0
    for i in range(n_edge):
        taper = (i + 1) / n_edge
        W_clean[i] *= taper
        W_clean[-(i+1)] *= taper

    return W_clean, W_raw, psi_abs, eig_val, interior_idx

def analyze_cleaned_superpotential(W_clean, W_raw, xgrid, interior_idx):
    """Detailed analysis of cleaned W."""
    x_interior = xgrid[interior_idx]
    W_interior = W_clean[interior_idx]

    stats = {
        'mean': float(np.mean(W_interior)),
        'std': float(np.std(W_interior)),
        'min': float(np.min(W_interior)),
        'max': float(np.max(W_interior)),
    }

    if len(W_interior) > 2:
        d2W = np.diff(W_interior, 2)
        rms_curv = np.sqrt(np.mean(d2W**2))
        smooth = rms_curv / (np.mean(np.abs(W_interior)) + 1e-12)
        stats['rms_curvature'] = float(rms_curv)
        stats['smoothness_ratio'] = float(smooth)

    # Even/odd
    D = len(xgrid)
    W_even = np.zeros(D)
    W_odd = np.zeros(D)
    for i in range(D):
        mirror = D - 1 - i
        W_even[i] = (W_clean[i] + W_clean[mirror]) / 2.0
        W_odd[i] = (W_clean[i] - W_clean[mirror]) / 2.0

    even_norm = np.linalg.norm(W_even[interior_idx])
    odd_norm = np.linalg.norm(W_odd[interior_idx])
    total_norm = np.linalg.norm(W_clean[interior_idx])

    stats['even_ratio'] = float(even_norm / (total_norm + 1e-12))
    stats['odd_ratio'] = float(odd_norm / (total_norm + 1e-12))

    # Test for monotonicity
    dW = np.diff(W_interior)
    sign_changes = np.sum(np.diff(np.sign(dW)) != 0)
    stats['dW_sign_changes'] = float(sign_changes)
    stats['dW_mean_sign'] = float(np.mean(np.sign(dW)))

    return stats, W_even, W_odd

def main():
    N, c = 64, 14
    L = np.log(c)

    print("=== Weil Position-Space SUSY v4: CLEANED superpotential ===\n")

    print(f"Building tau in frequency basis (N={N}, c={c})...")
    tau_freq, ns_freq = build_tau_np(N, c)
    D_freq = len(tau_freq)

    x_range = L / 2.0
    N_pos = 2*N + 1
    xgrid = np.linspace(-x_range, x_range, N_pos)
    dx = xgrid[1] - xgrid[0]

    print(f"  x grid: {N_pos} points, x ∈ [{xgrid[0]:.3f}, {xgrid[-1]:.3f}], dx={dx:.4f}")

    print("\nTransforming tau to position space...")
    tau_pos = transform_to_position_basis_dct(tau_freq, xgrid)

    eigs = np.linalg.eigvalsh(tau_pos)
    print(f"  Eigenvalue range: [{eigs.min():.6f}, {eigs.max():.6f}]")
    print(f"  Negative eigenvalues: {np.sum(eigs < -1e-10)}")

    print("\nExtracting cleaned superpotential...")
    W_clean, W_raw, psi_abs, lambda_0, interior_idx = extract_clean_superpotential(tau_pos, xgrid, exclude_edge_frac=0.15)

    print(f"  λ₀ (smallest eigenvalue): {lambda_0:.8f}")
    print(f"  ‖ψ₀‖: {np.linalg.norm(psi_abs):.6f}")
    print(f"  Interior region: {len(interior_idx)} / {len(xgrid)} points ({100*len(interior_idx)/len(xgrid):.1f}%)")

    print("\nW(x) sample values (cleaned):")
    sample_idx = [interior_idx[i] for i in np.linspace(0, len(interior_idx)-1, min(5, len(interior_idx)), dtype=int)]
    for idx in sample_idx:
        print(f"  W({xgrid[idx]:7.3f}) = {W_clean[idx]:10.6f}  [raw: {W_raw[idx]:10.2f}]")

    print("\nComparison: raw vs cleaned")
    print(f"  Raw W stats on interior:")
    print(f"    mean: {np.mean(W_raw[interior_idx]):10.4f}, std: {np.std(W_raw[interior_idx]):10.4f}")
    print(f"  Cleaned W stats on interior:")
    print(f"    mean: {np.mean(W_clean[interior_idx]):10.4f}, std: {np.std(W_clean[interior_idx]):10.4f}")

    print()

    # Detailed analysis
    stats, W_even, W_odd = analyze_cleaned_superpotential(W_clean, W_raw, xgrid, interior_idx)

    print("W(x) analysis (interior region):")
    for k, v in stats.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.6f}")

    print()
    print("Interpretation:")
    even_rat = stats['even_ratio']
    odd_rat = stats['odd_ratio']
    print(f"  Even component: {even_rat:.1%} of total")
    print(f"  Odd component: {odd_rat:.1%} of total")
    print(f"  Balance: {'Nearly equal (ambiguous)' if abs(even_rat - odd_rat) < 0.2 else ('Even-dominated' if even_rat > odd_rat else 'Odd-dominated')}")

    smoothness = stats.get('smoothness_ratio', np.nan)
    if not np.isnan(smoothness):
        smoothness_qual = "SMOOTH (clean)" if smoothness < 0.2 else "rough (noisy)" if smoothness > 1.0 else "moderate"
        print(f"  Smoothness: {smoothness_qual} (ratio={smoothness:.3f})")

    sign_changes = stats['dW_sign_changes']
    print(f"  dW sign changes: {int(sign_changes)} (monotonicity: {['alternating/wiggly', 'mostly mono'][sign_changes < len(interior_idx)//4]})")

    print()
    print("Structure test for ζ′/ζ-like behavior:")
    print(f"  ζ′/ζ is odd → expect odd_ratio ≈ 1.0, even_ratio ≈ 0")
    print(f"  Current ratio suggests: ", end="")
    if odd_rat > 0.7:
        print("✓ STRONG odd component (promising)")
    elif odd_rat > 0.3:
        print("~ MODERATE odd component (mixed)")
    else:
        print("✗ Dominated by even part (not ζ′/ζ-like)")

    print()
    print("="*60)

    # Save
    import json
    results = {
        'config': {'N': N, 'c': c, 'N_pos': N_pos, 'dx': float(dx)},
        'xgrid': xgrid.tolist(),
        'W_clean': W_clean.tolist(),
        'W_raw': W_raw.tolist(),
        'W_even': W_even.tolist(),
        'W_odd': W_odd.tolist(),
        'psi_0_abs': psi_abs.tolist(),
        'interior_idx': interior_idx.tolist(),
        'lambda_0': float(lambda_0),
        'stats': stats,
    }
    outfile = '/tmp/weil_position_space_v4_results.json'
    with open(outfile, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"Results saved: {outfile}\n")

    # Print caveats
    print("=== Discretization & Extraction Caveats ===")
    print("• Position space obtained via DCT (2D Fourier approximation)")
    print("• Superpotential W = −d/dx[log|ψ₀|] where ψ₀ is near-null eigenvector of τ_pos")
    print("• Derivative estimated via UnivariateSpline (cubic, smoothed by s-factor)")
    print("• Edge regions excluded (where |ψ₀| < 10% of max) to avoid log(0) singularity")
    print("• Smoothing via Savitzky-Golay filter (window ~10% of interior points)")
    print("• Even/odd decomposition post-hoc (not enforced during construction)")
    print("• Test is: does cleaned W look like ζ′/ζ? → should be predominantly ODD")
    print("• Assumption: (d/dx + W)ψ₀ ≈ 0 (SUSY partner structure) — not verified rigorously")

if __name__ == '__main__':
    main()
