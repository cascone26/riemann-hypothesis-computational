#!/usr/bin/env python3
"""ATTACK B3, first pass: Weil form in POSITION (log-line) space, SUSY superpotential extraction.

Construct a matrix on a log-space grid:
  - Kernel: Ω̌(x) = −e^{−|x|/2}/(1−e^{−2|x|})  [smooth, convolution on log-line]
  - Prime sector: point masses at x=±log(p) with weight Λ(p)p^{-1/2}
  - Boundary rank-2

Take matrix square root A (=A†, valid SUSY structure), extract superpotential W(x) from
near-diagonal (A ≈ d/dx + W(x) in discretized form).

TEST: Is W(x) clean/arithmetic in position space? Does it resemble −(log ψ₀)′, i.e., a ζ′/ζ-type object?

STATUS: First pass, partial acceptable. Report sample W values, smoothness, recognizability.
Report exact command and RAW output. Honesty on discretization limits.
"""
import numpy as np
from scipy.special import digamma
import sys
sys.path.insert(0, "operators")
from critical_eigvec_id import prime_powers_up_to

def von_mangoldt(n):
    """Von Mangoldt function Λ(n): log(p) if n=p^k, else 0."""
    if n <= 1:
        return 0.0
    # Factor n
    x = n
    for p in range(2, int(np.sqrt(n)) + 1):
        if x % p == 0:
            lam = np.log(p)
            while x % p == 0:
                x //= p
            if x == 1:  # n was a prime power p^k
                return lam
            else:  # composite, not a prime power
                return 0.0
    # n is prime
    return np.log(n)

def weil_kernel(x):
    """Ω̌(x) = −e^{−|x|/2} / (1 − e^{−2|x|}). Smooth on ℝ, even, arch-type."""
    ax = np.abs(x)
    if ax < 1e-12:
        return -0.25  # limit as x->0: -e^0/(1-e^0) ~ -1/(4) near the discontinuity
    num = -np.exp(-ax / 2.0)
    den = 1.0 - np.exp(-2.0 * ax)
    return num / den if den != 0 else 0.0

# vectorize
weil_kernel_vec = np.vectorize(weil_kernel)

def build_weil_position_matrix(xgrid, prime_powers, c_max):
    """
    Build Weil form matrix in position space.

    xgrid: log-space positions (np.array)
    prime_powers: list of (q,p) tuples from prime_powers_up_to(c_max)
    c_max: cutoff for primes

    Returns: symmetric matrix T[i,j] = convolution kernel + prime point masses
    """
    D = len(xgrid)
    T = np.zeros((D, D))
    dx = xgrid[1] - xgrid[0] if D > 1 else 1.0  # grid spacing

    # Arch kernel: convolution Ω̌(x_i - x_j)
    for i in range(D):
        for j in range(D):
            dx_ij = xgrid[i] - xgrid[j]
            T[i, j] += weil_kernel(dx_ij)

    # Boundary rank-2: add a mild symmetric term (finite-size effect)
    # T[0,0] += 1.0; T[-1,-1] += 1.0  # or leave as is for purity

    # Prime sector: point masses at ±log(q) weighted by Λ(p)p^{-1/2}
    for q, p in prime_powers:
        log_q = np.log(q)
        weight = np.log(p) * (q ** -0.5)

        # Find indices closest to ±log(q)
        for sign in [1, -1]:
            target_x = sign * log_q
            # Find nearest grid point
            idx = np.argmin(np.abs(xgrid - target_x))
            if 0 <= idx < D:
                # Add point mass contribution: diagonal + small off-diag structure
                T[idx, idx] += weight

    return T

def sqrtm_psd(M, eps=1e-12):
    """Compute matrix square root of PSD matrix via eigendecomposition."""
    w, V = np.linalg.eigh(M)
    w = np.maximum(w, eps)  # clip small negatives
    return V @ np.diag(np.sqrt(w)) @ V.T

def extract_superpotential_discrete(A, xgrid):
    """
    Extract superpotential W(x) from near-diagonal of A.

    For A ≈ d/dx + W(x) in discretized form:
    - Diagonal A[i,i] roughly encodes W(x_i) after removing kinetic part
    - Near-diagonal (±1) encodes the discrete derivative scaling

    Simple extraction: W(x_i) ≈ A[i,i] / dx, or fit from pattern
    """
    D = len(A)
    dx = xgrid[1] - xgrid[0] if D > 1 else 1.0

    # Read diagonal
    diag = np.diag(A)

    # Estimate derivative part from off-diag, remove it
    # A[i,i±1] ~ ±(1/2dx) from centered difference of d/dx
    # True W ≈ diag - (off-diag pattern effects)

    # Simplest: just report the diagonal itself (it contains both kinetic + potential)
    # For a true first-order operator d/dx + W(x):
    #   (d/dx + W)^2 = d^2/dx^2 + (dW/dx + W d/dx + W^2)
    # So the diagonal of the square involves W^2 and dW/dx terms

    # Better: fit W by noting A[i,i] should be the anti-Hermitian part + potential
    # For now, report diag directly
    return diag

def analyze_superpotential(W_values, xgrid):
    """Analyze smoothness and structure of superpotential."""
    # 2nd order finite difference (smoothness)
    if len(W_values) > 2:
        d2W = np.diff(W_values, 2)
        rms_curvature = np.sqrt(np.mean(d2W ** 2))
        mean_W = np.mean(np.abs(W_values[W_values != 0]))
        smoothness = rms_curvature / (mean_W + 1e-12)
    else:
        smoothness = np.nan

    # Check for monotonicity / structure
    dW = np.diff(W_values)
    sign_changes = np.sum(np.diff(np.sign(dW)) != 0)

    return {
        'mean': np.mean(W_values),
        'std': np.std(W_values),
        'rms_curvature': rms_curvature if len(W_values) > 2 else np.nan,
        'smoothness_ratio': smoothness,
        'sign_changes_in_dW': sign_changes,
        'min': np.min(W_values),
        'max': np.max(W_values),
    }

def main():
    # ===== Configuration =====
    c_max = 14  # prime cutoff, reasonable for first pass
    N_grid = 64  # points on log-line
    x_range = np.log(c_max)  # symmetric around 0

    xgrid = np.linspace(-x_range, x_range, N_grid)
    prime_powers = prime_powers_up_to(c_max)

    print(f"=== Weil Position-Space SUSY Superpotential (c_max={c_max}, N={N_grid}) ===\n")
    print(f"Log-space grid: x ∈ [{xgrid[0]:.3f}, {xgrid[-1]:.3f}], dx ≈ {xgrid[1]-xgrid[0]:.4f}")
    print(f"Prime powers up to {c_max}: {len(prime_powers)} total")
    print(f"  First 10: {prime_powers[:10]}")
    print()

    # ===== Build Weil matrix =====
    print("Building Weil form matrix (arch kernel + prime point masses)...")
    T = build_weil_position_matrix(xgrid, prime_powers, c_max)

    # Check symmetry
    sym_err = np.linalg.norm(T - T.T) / (np.linalg.norm(T) + 1e-12)
    print(f"  Symmetry error: {sym_err:.2e}")

    # Check PSD (for valid matrix square root)
    eigs_T = np.linalg.eigvalsh(T)
    n_neg = np.sum(eigs_T < -1e-10)
    print(f"  Eigenvalue range: [{eigs_T.min():.4f}, {eigs_T.max():.4f}]")
    print(f"  Negative eigenvalues: {n_neg} (if >0, matrix is not PSD; square root may be complex)")
    print()

    # ===== Matrix square root =====
    print("Computing matrix square root A = √T...")
    A = sqrtm_psd(T)

    # Verify A^2 ≈ T
    A2 = A @ A
    err_sq = np.linalg.norm(A2 - T) / (np.linalg.norm(T) + 1e-12)
    print(f"  Reconstruction error ‖A² - T‖ / ‖T‖ = {err_sq:.2e}")

    # Check symmetry of A
    sym_A = np.linalg.norm(A - A.T) / (np.linalg.norm(A) + 1e-12)
    print(f"  A symmetry error: {sym_A:.2e} (should be ~0 for real-symmetric A)")
    print()

    # ===== Extract superpotential =====
    print("Extracting superpotential W(x) from near-diagonal of A...")
    W = extract_superpotential_discrete(A, xgrid)

    # Report sample values
    print(f"\nSuperpotential W(x) sample values:")
    print(f"  W(x_{{center}}): {W[N_grid//2]:.6f}")
    sample_indices = [0, N_grid//4, N_grid//2, 3*N_grid//4, N_grid-1]
    for idx in sample_indices:
        if 0 <= idx < N_grid:
            print(f"    W(x[{idx}]={xgrid[idx]:7.3f}) = {W[idx]:10.6f}")

    # Analyze smoothness
    stats = analyze_superpotential(W, xgrid)
    print(f"\nSuperpotential statistics:")
    for key, val in stats.items():
        if np.isnan(val):
            print(f"  {key}: NaN")
        else:
            print(f"  {key}: {val:.6f}")

    print()

    # ===== Check for arithmetic structure =====
    print("Testing for ζ′/ζ-type structure:")
    # A ζ′/ζ function is logarithmic derivative of a Dirichlet series
    # Its real part is often non-monotonic but "smooth on average"
    # Check: is W close to a derivative-of-log function?

    # Try fitting W ≈ d/dx[log ψ(x)] for some ψ
    # This is hard without knowing ψ, but we can check:
    # - Is W anti-symmetric? (ζ′/ζ is odd)
    # - Is |W| monotonic or wiggly?
    # - Is W close to a rational function or transcendental?

    W_center = W[N_grid//2 : N_grid//2+10]
    W_antisym_test = np.array([W[i] + W[N_grid-1-i] for i in range(N_grid//2)])
    antisym_deviation = np.std(W_antisym_test) / (np.std(W) + 1e-12)

    print(f"  Anti-symmetry test (W(-x) + W(x) should be small): {antisym_deviation:.4f}")
    print(f"  Interpretation: {'weak' if antisym_deviation < 0.2 else 'strong'} anti-symmetric component")

    # Derivative-like structure: is dW/dx well-defined?
    if len(W) > 2:
        dW = np.diff(W)
        d2W = np.diff(dW)
        print(f"  dW/dx roughness (RMS of d²W): {np.sqrt(np.mean(d2W**2)):.6f}")
        print(f"  dW/dx mean sign bias: {np.mean(np.sign(dW)):.3f} (0=neutral, ±1=mono)")

    print()

    # ===== Discretization caveats =====
    print("=== Discretization Caveats & Honesty ===")
    print(f"• Grid spacing dx = {xgrid[1]-xgrid[0]:.4f} (position space): finite-difference")
    print(f"  derivative estimates A[i,i±1] smear the true d/dx structure")
    print(f"• Weil kernel Ω̌(x) integrated only via point-evaluation on grid,")
    print(f"  not via full convolution over continuum")
    print(f"• Prime point masses added only to nearest grid point (quantization noise)")
    print(f"• Matrix square root via eigendecomposition: stable numerically, but")
    print(f"  suppresses high-frequency (high-eigenvalue) artifacts")
    print(f"• Superpotential extraction is NAIVE: A[i,i] without deconvolution of kinetic part")
    print(f"  → true W(x) requires subtracting d/dx contribution, only done implicitly")
    print()

    # ===== Save results =====
    results = {
        'xgrid': xgrid.tolist(),
        'W': W.tolist(),
        'A': A.tolist(),
        'T': T.tolist(),
        'stats': {k: (float(v) if not np.isnan(v) else None) for k,v in stats.items()},
        'antisym_deviation': float(antisym_deviation),
        'config': {
            'c_max': c_max,
            'N_grid': N_grid,
            'x_range': float(x_range),
        }
    }

    import json
    outfile = "/tmp/weil_position_space_susy_results.json"
    with open(outfile, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Results saved to {outfile}")

if __name__ == '__main__':
    main()
