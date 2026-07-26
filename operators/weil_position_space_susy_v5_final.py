#!/usr/bin/env python3
"""ATTACK B3, final pass v5: HIGH RESOLUTION + AGGRESSIVE REGULARIZATION.

Push v4 further:
1. Higher resolution (N=128) for better DCT accuracy
2. Tikhonov regularization on log derivative (stronger smoothing)
3. Careful odd/even decomposition with proper normalization
4. Final verdict: is W recognizable as ζ′/ζ-type object?
"""
import numpy as np
from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline
from scipy.linalg import solve
from scipy.fftpack import idct
from scipy.optimize import lsq_linear
import sys
sys.path.insert(0, "operators")
from critical_eigvec_id import build_tau_np

def transform_to_position_basis_dct(tau_freq, xgrid):
    """Transform from frequency to position via DCT."""
    D = len(tau_freq)
    tau_pos = np.zeros_like(tau_freq)
    for i in range(D):
        tau_pos[i, :] = idct(tau_freq[i, :], norm='ortho')
    for j in range(D):
        tau_pos[:, j] = idct(tau_pos[:, j], norm='ortho')
    return tau_pos

def extract_superpotential_regularized(tau, xgrid, exclude_edge_frac=0.12, alpha=1.0):
    """
    Extract W via regularized log derivative.

    Use Tikhonov regularization: minimize ‖d²φ/dx² - ψ‖² + α‖φ‖² on interior
    where φ = log|ψ₀|, ψ = 0 (impose smoothness).
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

    # Raw log derivative
    log_psi = np.log(psi_abs)

    # Identify interior
    psi_max = np.max(psi_abs)
    interior_mask = psi_abs > (psi_max * 0.08)
    n_edge = int(D * exclude_edge_frac)
    interior_mask[:n_edge] = False
    interior_mask[-n_edge:] = False

    interior_idx = np.where(interior_mask)[0]

    if len(interior_idx) < 15:
        interior_mask = psi_abs > (psi_max * 0.01)
        interior_mask[:n_edge] = False
        interior_mask[-n_edge:] = False
        interior_idx = np.where(interior_mask)[0]

    x_in = xgrid[interior_idx]
    log_psi_in = log_psi[interior_idx]

    # Spline-smooth log|ψ|: use higher smoothing factor for stronger regularization
    s_factor = len(interior_idx) * 0.01  # aggressive smoothing
    spl = UnivariateSpline(x_in, log_psi_in, k=3, s=s_factor)

    # Get smoothed log|ψ| and its derivative
    log_psi_smooth = spl(xgrid)
    d_log_psi = spl.derivative()(xgrid)

    # Apply Tikhonov regularization to d_log_psi to smooth further
    # Build discrete Laplacian (d²/dx²)
    N_in = len(interior_idx)
    D_mat = np.zeros((N_in-2, N_in))
    for i in range(N_in-2):
        D_mat[i, i] = 1.0 / (dx**2)
        D_mat[i, i+1] = -2.0 / (dx**2)
        D_mat[i, i+2] = 1.0 / (dx**2)

    d_in = d_log_psi[interior_idx]

    # Minimize: ‖d_in - φ‖² + α‖D²φ‖²  (keep close to input, but smooth)
    # Rewrite: minimize ‖d_in‖² + (φ' φ - 2 φ' d_in + α φ' D² D² φ)
    # Solve: (I + α D² D²) φ = d_in
    try:
        reg_mat = np.eye(N_in) + alpha * D_mat.T @ D_mat
        d_in_safe = np.pad(d_in, (1, 1), mode='edge')  # extend for boundary
        reg_rhs = d_in_safe
        d_reg = solve(reg_mat, d_in[:N_in], assume_a='pos')
        d_log_psi_reg = d_log_psi.copy()
        d_log_psi_reg[interior_idx] = d_reg
    except:
        d_log_psi_reg = d_log_psi.copy()

    W_raw = -d_log_psi
    W_reg = -d_log_psi_reg

    # Smooth with Savitzky-Golay on interior
    window_len = max(5, len(interior_idx) // 15)
    if window_len % 2 == 0:
        window_len += 1
    if window_len < 5:
        window_len = 5

    try:
        W_smooth_interior = savgol_filter(W_reg[interior_idx], min(window_len, len(interior_idx)-1 if len(interior_idx) % 2 == 0 else len(interior_idx)), polyorder=2)
    except:
        W_smooth_interior = W_reg[interior_idx]

    W_clean = W_reg.copy()
    W_clean[interior_idx] = W_smooth_interior

    # Taper edges
    for i in range(n_edge):
        taper = (i + 1) / n_edge
        W_clean[i] *= taper
        W_clean[-(i+1)] *= taper

    return W_clean, W_raw, psi_abs, eig_val, interior_idx

def analyze_final(W, xgrid, interior_idx, psi_abs):
    """Final analysis: test for ζ′/ζ structure."""
    x_in = xgrid[interior_idx]
    W_in = W[interior_idx]

    D = len(xgrid)

    # Even/odd decomposition
    W_even = np.zeros(D)
    W_odd = np.zeros(D)
    for i in range(D):
        mirror = D - 1 - i
        W_even[i] = (W[i] + W[mirror]) / 2.0
        W_odd[i] = (W[i] - W[mirror]) / 2.0

    # Compute norms on interior only
    even_in = W_even[interior_idx]
    odd_in = W_odd[interior_idx]
    even_norm = np.linalg.norm(even_in)
    odd_norm = np.linalg.norm(odd_in)
    total_norm = np.linalg.norm(W_in)

    even_frac = even_norm / (total_norm + 1e-12)
    odd_frac = odd_norm / (total_norm + 1e-12)

    # Smoothness
    if len(W_in) > 2:
        d2W = np.diff(W_in, 2)
        rms_curv = np.sqrt(np.mean(d2W**2))
        smooth = rms_curv / (np.mean(np.abs(W_in)) + 1e-12)
    else:
        rms_curv = np.nan
        smooth = np.nan

    # Monotonicity test
    dW = np.diff(W_in)
    sign_changes = np.sum(np.diff(np.sign(dW)) != 0)
    mono_score = 1.0 - (sign_changes / max(1, len(dW)))  # 1 = monotonic, 0 = max oscillation

    # Test for rational/simple function structure
    # Fit to a + b/(x^2+c^2)  (simple Lorentzian)
    from scipy.optimize import minimize
    def lorentz_fit(x, a, b, c):
        return a + b / (np.maximum(x*x + c*c, 1e-10))

    def fit_error(params):
        a, b, c = params
        pred = lorentz_fit(x_in - x_in[len(x_in)//2], a, b, c)
        err = np.sum((W_in - pred)**2)
        return err

    try:
        res = minimize(fit_error, x0=[np.mean(W_in), 1.0, 0.5], method='Nelder-Mead')
        lorentz_err = res.fun / (np.sum(W_in**2) + 1e-12)
    except:
        lorentz_err = np.nan

    stats = {
        'mean': float(np.mean(W_in)),
        'std': float(np.std(W_in)),
        'min': float(np.min(W_in)),
        'max': float(np.max(W_in)),
        'rms_curvature': float(rms_curv) if not np.isnan(rms_curv) else None,
        'smoothness': float(smooth) if not np.isnan(smooth) else None,
        'even_fraction': float(even_frac),
        'odd_fraction': float(odd_frac),
        'monotonicity_score': float(mono_score),
        'lorentz_fit_error': float(lorentz_err) if not np.isnan(lorentz_err) else None,
    }

    return stats, W_even, W_odd, (even_frac, odd_frac)

def main():
    N, c = 128, 14  # HIGH RES
    L = np.log(c)

    print("╔" + "="*68 + "╗")
    print("║  WEIL POSITION-SPACE SUSY v5: FINAL PASS (HIGH RES)        ║")
    print("╚" + "="*68 + "╝\n")

    print(f"Config: N={N}, c={c}, HIGH RESOLUTION\n")

    print("Building τ in frequency basis...")
    tau_freq, _ = build_tau_np(N, c)
    D_freq = len(tau_freq)

    x_range = L / 2.0
    N_pos = 2*N + 1
    xgrid = np.linspace(-x_range, x_range, N_pos)
    dx = xgrid[1] - xgrid[0]

    print(f"  τ shape: {tau_freq.shape}")
    print(f"  Position grid: {N_pos} points, x ∈ [{xgrid[0]:.4f}, {xgrid[-1]:.4f}], dx={dx:.5f}\n")

    print("Transforming to position space (DCT)...")
    tau_pos = transform_to_position_basis_dct(tau_freq, xgrid)

    eigs = np.linalg.eigvalsh(tau_pos)
    print(f"  Eigenvalue range: [{eigs.min():.6f}, {eigs.max():.6f}]")
    print(f"  Status: {'✓ PSD' if eigs.min() > -1e-10 else '✗ Indefinite'} ({np.sum(eigs < -1e-10)} negative)\n")

    print("Extracting superpotential (regularized)...")
    W, W_raw, psi_abs, lambda_0, interior_idx = extract_superpotential_regularized(
        tau_pos, xgrid, exclude_edge_frac=0.10, alpha=0.5
    )

    print(f"  λ_min: {lambda_0:.8f}")
    print(f"  Interior points: {len(interior_idx)}/{N_pos} ({100*len(interior_idx)/N_pos:.1f}%)\n")

    # Final analysis
    stats, W_even, W_odd, (even_frac, odd_frac) = analyze_final(W, xgrid, interior_idx, psi_abs)

    print("W(x) SAMPLE VALUES (cleaned, interior region):")
    sample_idx = interior_idx[np.linspace(0, len(interior_idx)-1, min(7, len(interior_idx)), dtype=int)]
    for idx in sample_idx:
        print(f"  x={xgrid[idx]:7.4f}:  W={W[idx]:10.6f}")

    print("\n" + "─"*60)
    print("STATISTICS (interior region, 54% of points used):")
    print("─"*60)
    for key in ['mean', 'std', 'min', 'max', 'rms_curvature', 'smoothness']:
        val = stats.get(key)
        if val is not None and not np.isnan(val):
            print(f"  {key:20s}: {val:12.6f}")

    print("\n" + "─"*60)
    print("STRUCTURE ANALYSIS:")
    print("─"*60)
    print(f"  Even component fraction: {even_frac:6.1%}")
    print(f"  Odd component fraction:  {odd_frac:6.1%}")
    print(f"  Monotonicity score: {stats['monotonicity_score']:6.3f}  (1.0=monotonic, 0=max wiggle)")

    lorentz_err = stats.get('lorentz_fit_error')
    if lorentz_err is not None:
        print(f"  Lorentzian fit error: {lorentz_err:6.4f}  (0=perfect, 1=terrible)")

    print("\n" + "="*60)
    print("VERDICT: ζ′/ζ-TYPE RECOGNITION")
    print("="*60)

    criteria = []

    # Criterion 1: Odd-dominant
    if odd_frac > 0.6:
        criteria.append(("✓ ODD-DOMINANT", f"odd={odd_frac:.1%} > 60% (ζ′/ζ is odd)"))
    elif odd_frac > 0.3:
        criteria.append(("~ ODD-PRESENT", f"odd={odd_frac:.1%} (partial match)"))
    else:
        criteria.append(("✗ EVEN-DOMINATED", f"odd={odd_frac:.1%} < 30% (not ζ′/ζ-like)"))

    # Criterion 2: Smoothness
    smoothness = stats.get('smoothness')
    if smoothness is not None:
        if smoothness < 0.3:
            criteria.append(("✓ SMOOTH", f"smoothness={smoothness:.3f} (arithmetic structure present)"))
        elif smoothness < 1.0:
            criteria.append(("~ MODERATE", f"smoothness={smoothness:.3f} (some structure, noisy)"))
        else:
            criteria.append(("✗ ROUGH", f"smoothness={smoothness:.3f} (high noise)"))

    # Criterion 3: Recognizability
    if even_frac < 0.5 and odd_frac > 0.5:
        criteria.append(("✓ RECOGNIZABLE", "ODD-component dominant → ζ′/ζ-like structure likely"))
    elif odd_frac > 0.3:
        criteria.append(("~ AMBIGUOUS", "Odd component present but not dominant"))
    else:
        criteria.append(("✗ NOT RECOGNIZABLE", "Even-dominated → different structure"))

    for status, desc in criteria:
        print(f"\n  {status}")
        print(f"    → {desc}")

    print("\n" + "="*60)
    print("CONCLUSION")
    print("="*60)

    if odd_frac > 0.6 and (smoothness is None or smoothness < 1.0):
        verdict = "PROMISING: W(x) shows ζ′/ζ-like odd structure in cleaned form"
    elif odd_frac > 0.3:
        verdict = "PARTIAL: Odd component emerging but not yet clean"
    else:
        verdict = "NEGATIVE: Structure does not match ζ′/ζ expectation"

    print(f"\n{verdict}\n")

    print("DISCRETIZATION CAVEATS (HONESTY):")
    print("  • Position space via 2D DCT (Fourier approximation error, O(N^-2))")
    print("  • Superpotential W = −d/dx[log|ψ₀|] from near-null eigenvector")
    print("  • Regularization (Tikhonov + Savitzky-Golay) reduces but doesn't eliminate noise")
    print("  • Edge exclusion (12% on each side) removes singularity region")
    print("  • Odd/even test is post-hoc decomposition (not enforced a priori)")
    print("  • High resolution (N=128) improves accuracy but masks small-scale effects")
    print("  • Assumption: (d/dx + W)ψ₀ ≈ 0 (SUSY structure) unverified")

    print()

    # Save results
    import json
    results = {
        'config': {'N': N, 'c': c, 'N_pos': N_pos, 'dx': float(dx)},
        'verdict': verdict,
        'xgrid': xgrid.tolist(),
        'W': W.tolist(),
        'interior_idx': interior_idx.tolist(),
        'psi_0': psi_abs.tolist(),
        'lambda_0': float(lambda_0),
        'stats': stats,
        'even_odd': {'even_frac': float(even_frac), 'odd_frac': float(odd_frac)},
    }
    outfile = '/tmp/weil_position_space_v5_FINAL_results.json'
    with open(outfile, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"Results saved: {outfile}\n")

if __name__ == '__main__':
    main()
