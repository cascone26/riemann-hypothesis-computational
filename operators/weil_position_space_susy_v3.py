#!/usr/bin/env python3
"""ATTACK B3, third pass: transform EXISTING tau to position space via Fourier.

Instead of reconstructing the Weil form naively, leverage the KNOWN tau from
critical_eigvec_id.build_tau_np (frequency basis), then change basis to position
space via DCT or explicit Fourier. The frequency-basis tau is already tuned
(arch+boundary+primes integrated properly in frequency domain).

Hypothesis: position-space view of the same tau may reveal cleaner arithmetic
structure in W(x).
"""
import numpy as np
from scipy.special import digamma
from scipy.fftpack import dct, idct
import sys
sys.path.insert(0, "operators")
from critical_eigvec_id import build_tau_np

def transform_to_position_basis_dct(tau_freq, xgrid):
    """
    Approximate position-space view via Discrete Cosine Transform.

    tau_freq is (2N+1, 2N+1) symmetric in Fourier basis.
    xgrid is position-space grid (assumed symmetric around 0).

    Use DCT to approximate: position-space view ≈ DCT_inverse(tau_freq).
    """
    D = len(tau_freq)
    # tau_freq is in index basis (-N:N+1), symmetric
    # Apply DCT-II to each row (approximates Fourier to position)
    tau_pos = np.zeros_like(tau_freq)
    for i in range(D):
        tau_pos[i, :] = idct(tau_freq[i, :], norm='ortho')
    # Also columns
    for j in range(D):
        tau_pos[:, j] = idct(tau_pos[:, j], norm='ortho')

    return tau_pos

def sqrtm_svd_clipped(M, k_clip=None):
    """
    SVD-based square root, keeping top k components.

    For indefinite M with large dynamic range, keep only largest singular
    values to avoid numeric issues.
    """
    U, s, Vt = np.linalg.svd(M)
    if k_clip is None:
        # Keep 95% of singular value energy
        s_energy = np.cumsum(s**2) / np.sum(s**2)
        k_clip = np.argmax(s_energy >= 0.95) + 1
    print(f"  SVD truncation: keeping {k_clip} / {len(s)} singular values")
    s_clip = s[:k_clip]
    U_clip = U[:, :k_clip]
    Vt_clip = Vt[:k_clip, :]

    s_sqrt = np.sqrt(np.maximum(s_clip, 1e-14))
    A_approx = U_clip @ np.diag(s_sqrt) @ Vt_clip
    return A_approx

def extract_superpotential_from_eigenvector(tau, xgrid):
    """
    Extract W(x) from the NEAR-NULL eigenvector of tau.

    For a SUSY operator H = d²/dx² + dW/dx + W², the ground state ψ₀
    satisfies: (d/dx + W) ψ₀ = 0  =>  W = −(d/dx log ψ₀)′  wait, wrong sign
    Actually: (d/dx + W) ψ₀ = 0  =>  W = −ψ₀'/ψ₀

    If ψ₀ is near-null eigenvector of tau, extract W ≈ −d/dx[log ψ₀].
    """
    w, Q = np.linalg.eigh(tau)
    idx_min = np.argmin(np.abs(w))
    psi = Q[:, idx_min]
    eig_val = w[idx_min]

    # Normalize
    psi = psi / np.linalg.norm(psi)

    D = len(xgrid)
    dx = xgrid[1] - xgrid[0] if D > 1 else 1.0

    # Compute d/dx[log|psi|] ≈ (psi'/psi)
    # Use finite differences (or complex log derivative)
    psi_abs = np.abs(psi)
    psi_abs = np.maximum(psi_abs, 1e-14)  # avoid log(0)

    # d/dx[log psi_abs] via finite diff
    log_psi = np.log(psi_abs)
    d_log_psi = np.diff(log_psi) / dx

    # Pad to match original length
    d_log_psi_full = np.pad(d_log_psi, (1, 0), mode='edge')  # forward diff at start

    # W ≈ −d²/dx²[log psi] / 2  or just  W ≈ −d/dx[log psi] depending on normalization
    # For simplicity: W ≈ −d/dx[log psi]
    W_from_psi = -d_log_psi_full

    return W_from_psi, psi, eig_val

def analyze_arithmetic(W, xgrid):
    """Test if W has clean arithmetic/transcendental structure."""
    D = len(W)
    dx = xgrid[1] - xgrid[0]

    # Even/odd decomposition
    W_even, W_odd = np.zeros(D), np.zeros(D)
    for i in range(D):
        mirror = D - 1 - i
        W_even[i] = (W[i] + W[mirror]) / 2.0
        W_odd[i] = (W[i] - W[mirror]) / 2.0

    even_ratio = np.linalg.norm(W_even) / (np.linalg.norm(W) + 1e-12)
    odd_ratio = np.linalg.norm(W_odd) / (np.linalg.norm(W) + 1e-12)

    # Test: fit to simple functions
    x_c = xgrid[D//2:]
    W_c = W[D//2:]

    # Fit to Lorentzian 1/(x^2 + a^2)
    def lorentz(x, a):
        return 1.0 / (np.maximum(x*x + a*a, 1e-10))

    from scipy.optimize import minimize_scalar
    def fit_error(a):
        if a < 0:
            return 1e10
        pred = lorentz(np.abs(x_c), a)
        pred = pred / np.max(pred)  # normalize
        W_norm = W_c / (np.max(np.abs(W_c)) + 1e-12)
        return np.sum((W_norm - pred) ** 2)

    res = minimize_scalar(fit_error, bounds=(0.1, 10), method='bounded')
    lorentz_fit_error = res.fun

    return {
        'even_ratio': float(even_ratio),
        'odd_ratio': float(odd_ratio),
        'lorentz_fit_error': float(lorentz_fit_error),
        'W_even': W_even,
        'W_odd': W_odd,
    }

def main():
    N, c = 64, 14  # freq-basis params
    L = np.log(c)

    print("=== Weil Position-Space SUSY v3: Transform from frequency basis ===\n")
    print(f"Building tau in frequency basis (N={N}, c={c})...")
    tau_freq, ns_freq = build_tau_np(N, c)
    D_freq = len(tau_freq)
    print(f"  tau shape: {tau_freq.shape}")

    # Position space grid (symmetric around 0)
    x_range = L / 2.0  # cover roughly the same range
    N_pos = 2*N + 1
    xgrid = np.linspace(-x_range, x_range, N_pos)
    dx = xgrid[1] - xgrid[0]

    print(f"  x grid: {N_pos} points, x ∈ [{xgrid[0]:.3f}, {xgrid[-1]:.3f}], dx={dx:.4f}\n")

    # Transform to position basis via DCT
    print("Transforming tau to position space via DCT...")
    tau_pos = transform_to_position_basis_dct(tau_freq, xgrid)
    print(f"  tau_pos shape: {tau_pos.shape}")

    sym_err = np.linalg.norm(tau_pos - tau_pos.T) / (np.linalg.norm(tau_pos) + 1e-12)
    print(f"  Symmetry error: {sym_err:.2e}")

    eigs = np.linalg.eigvalsh(tau_pos)
    n_neg = np.sum(eigs < -1e-10)
    print(f"  Eigenvalue range: [{eigs.min():.4f}, {eigs.max():.4f}]")
    print(f"  Negative eigenvalues: {n_neg}")
    print()

    # Extract superpotential from near-null eigenvector
    print("Extracting superpotential from near-null eigenvector of tau_pos...")
    W, psi_0, lambda_0 = extract_superpotential_from_eigenvector(tau_pos, xgrid)

    print(f"  Smallest eigenvalue of tau_pos: {lambda_0:.6f}")
    print(f"  ‖ψ₀‖: {np.linalg.norm(psi_0):.6f}")
    print(f"  |ψ₀| range: [{np.min(np.abs(psi_0)):.4e}, {np.max(np.abs(psi_0)):.4e}]")
    print()

    print(f"W(x) = −d/dx[log|ψ₀|] sample values:")
    for idx in [0, N_pos//4, N_pos//2, 3*N_pos//4, N_pos-1]:
        if idx < len(W):
            print(f"  W({xgrid[idx]:7.3f}) = {W[idx]:10.6f}")

    print()

    # Statistics
    print("W(x) statistics:")
    print(f"  mean: {np.mean(W):.6f}")
    print(f"  std: {np.std(W):.6f}")
    print(f"  min: {np.min(W):.6f}")
    print(f"  max: {np.max(W):.6f}")

    if len(W) > 2:
        d2W = np.diff(W, 2)
        rms_curv = np.sqrt(np.mean(d2W**2))
        smooth = rms_curv / (np.mean(np.abs(W)) + 1e-12)
        print(f"  RMS curvature (d²W): {rms_curv:.6f}")
        print(f"  Smoothness ratio: {smooth:.6f}")

    print()

    # Arithmetic analysis
    print("Arithmetic structure test:")
    arith = analyze_arithmetic(W, xgrid)
    print(f"  Even component ratio: {arith['even_ratio']:.4f}")
    print(f"  Odd component ratio: {arith['odd_ratio']:.4f}")
    print(f"  Lorentzian fit error: {arith['lorentz_fit_error']:.6f}")

    print()

    # Try to factorize tau_pos = A^T A via truncated SVD
    print("Testing SUSY factorization (truncated SVD)...")
    A_svd = sqrtm_svd_clipped(tau_pos, k_clip=int(D_freq * 0.7))
    A_squared = A_svd @ A_svd.T
    recon_err = np.linalg.norm(A_squared - tau_pos) / (np.linalg.norm(tau_pos) + 1e-12)
    print(f"  Reconstruction error ‖A²−τ‖/‖τ‖: {recon_err:.4f}")

    print()
    print("=== Caveats ===")
    print(f"• Position space obtained via DCT (Fourier approximation)")
    print(f"• Near-null eigenvector of tau_pos provides W via logarithmic derivative")
    print(f"• W is estimated via finite differences (d/dx[log|ψ₀|])")
    print(f"• tau_pos may not be fully PSD (indefinite form from Riemann structure)")
    print(f"• Superpotential extraction assumes (d/dx + W)ψ₀ ≈ 0 (not verified)")

    print()
    print("="*60)

    # Save
    import json
    results = {
        'config': {'N': N, 'c': c, 'N_pos': N_pos, 'dx': float(dx)},
        'xgrid': xgrid.tolist(),
        'W': W.tolist(),
        'psi_0_abs': np.abs(psi_0).tolist(),
        'lambda_0': float(lambda_0),
        'arith': arith,
        'stats': {
            'mean': float(np.mean(W)),
            'std': float(np.std(W)),
            'min': float(np.min(W)),
            'max': float(np.max(W)),
        }
    }
    outfile = '/tmp/weil_position_space_v3_results.json'
    with open(outfile, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved: {outfile}")

if __name__ == '__main__':
    main()
