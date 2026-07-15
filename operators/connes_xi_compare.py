"""
Direct comparison: ξ̂(z) vs Ξ(z) (Riemann Xi function)

Connes' claim: as λ→∞, the Fourier transform of the minimal eigenvector
of QW_λ converges to the Riemann Xi function Ξ(z).

For λ²=14 (primes ≤ 13), the paper claims the first zero is predicted
to 10^{-55} error. This should be visible as ξ̂(z)/Ξ(z) ≈ constant.

We test this by:
1. Building the Weil matrix using v3's correct archimedean formula
2. Finding min eigenvector ξ
3. Computing ξ̂(z) for z on the real line
4. Computing Ξ(z) via mpmath
5. Checking if the RATIO ξ̂(z)/Ξ(z) is approximately constant
6. Checking zeros of ξ̂ vs zeros of Ξ (= zeta zeros on critical line)
7. Testing for multiple λ values to see convergence rate
"""

import mpmath
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import json
import os
import time

mpmath.mp.dps = 50  # 50 decimal places

DATA_DIR    = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")

# ── Riemann Xi function ────────────────────────────────────────────────────
def Xi(z):
    """
    Ξ(z) = ξ(1/2 + iz)  where ξ(s) = (1/2)s(s-1)π^{-s/2}Γ(s/2)ζ(s)
    Real-valued on the real line. All zeros are on the real line iff RH.
    """
    s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(z))
    return mpmath.mpf('0.5') * s * (s-1) * mpmath.power(mpmath.pi, -s/2) * mpmath.gamma(s/2) * mpmath.zeta(s)


# ── Archimedean kernel functions (from v3) ──────────────────────────────────

def rho(x):
    """ρ(x) = e^{x/2} / (2sinh(x)), real x > 0"""
    if x < 1e-10:
        return 1.0 / (2 * x + 1e-300)
    return np.exp(x/2) / (np.exp(x) - np.exp(-x))


def alpha_L(n, L, epsabs=1e-12):
    """α_L(n) = (1/π) ∫₀ᴸ sin(2πnx/L) ρ(x) dx"""
    if n == 0:
        return 0.0
    a = 2 * np.pi * n / L
    def integrand(x):
        if x < 1e-12:
            return a / 2.0
        return np.sin(a * x) * rho(x)
    result, _ = quad(integrand, 0, L, limit=300, epsabs=epsabs, epsrel=epsabs,
                     points=[1e-8, 0.01, 0.1, L/2])
    return result / np.pi


def beta_L(n, L, epsabs=1e-12):
    """β_L(n) = (1/L) ∫₀ᴸ x cos(2πnx/L) ρ(x) dx"""
    a = 2 * np.pi * n / L
    def integrand(x):
        if x < 1e-12:
            return 0.5
        return x * np.cos(a * x) * rho(x)
    result, _ = quad(integrand, 0, L, limit=300, epsabs=epsabs, epsrel=epsabs,
                     points=[1e-8, 0.01, 0.1, L/2])
    return result / L


def gamma_L(n, L, epsabs=1e-12):
    """γ_L(n) = ∫₀ᴸ [cos(2πnx/L) - e^{-x/2}] ρ(x) dx  (without constant)"""
    a = 2 * np.pi * n / L
    def integrand(x):
        if x < 1e-12:
            return 1/4.0 - a**2/4.0 * x if x > 1e-300 else 1/4.0
        return (np.cos(a * x) - np.exp(-x/2)) * rho(x)
    result, _ = quad(integrand, 0, L, limit=300, epsabs=epsabs, epsrel=epsabs,
                     points=[1e-8, 0.01, 0.1, L/2])
    return result


def W_arch(n, m, L, c_L, alpha_cache=None):
    """Correct archimedean Weil term W_R(V_n, V_m)"""
    if n == m:
        g = gamma_L(n, L) + c_L
        b = beta_L(n, L)
        return 2*g - 2*b
    else:
        if alpha_cache is not None:
            a_m = alpha_cache.setdefault(m, alpha_L(m, L))
            a_n = alpha_cache.setdefault(n, alpha_L(n, L))
        else:
            a_m = alpha_L(m, L)
            a_n = alpha_L(n, L)
        return (a_m - a_n) / (n - m)


def W_boundary(n, m, L):
    L_mp = mpmath.mpf(str(L))
    pi   = float(mpmath.pi)
    num  = 32*L_mp * mpmath.sinh(L_mp/4)**2 * (L_mp**2 - 16*pi**2*m*n)
    den  = (L_mp**2 + 16*pi**2*m**2) * (L_mp**2 + 16*pi**2*n**2)
    if abs(float(den)) < 1e-50:
        return 0.0
    return float(num / den)


def W_prime(n, m, L, primes):
    pi    = float(mpmath.pi)
    k_max = int(float(mpmath.exp(mpmath.mpf(str(L))))) + 1
    total = 0.0
    for p in primes:
        pk = p
        while pk <= k_max:
            y      = np.log(float(pk))
            lam_k  = np.log(float(p))
            weight = lam_k * pk**(-0.5)
            if n != m:
                q = (np.sin(2*pi*m*y/L) - np.sin(2*pi*n*y/L)) / (pi*(n-m))
            else:
                q = 2*(1 - abs(y)/L) * np.cos(2*pi*n*y/L)
            total += weight * q
            pk    *= p
    return total


def build_weil_matrix(N, L, primes):
    """Build the full Weil quadratic form matrix."""
    dim = 2*N + 1
    tau = np.zeros((dim, dim))
    alpha_cache = {}

    # Calibrate c_L: W_R(V_0,V_0) = -γ_E - log(2)
    expected_0 = float(-mpmath.euler - mpmath.log(2))
    b0 = beta_L(0, L)
    g0 = gamma_L(0, L)
    c_L = (expected_0 - 2*g0 + 2*b0) / 2.0

    t0 = time.time()
    print(f"  Building {dim}×{dim} Weil matrix (c_L={c_L:.6f})...", flush=True)
    for i in range(dim):
        n = i - N
        for j in range(i, dim):
            m = j - N
            arch = W_arch(n, m, L, c_L, alpha_cache)
            bnd  = W_boundary(n, m, L)
            prm  = W_prime(n, m, L, primes)
            val  = arch + bnd - prm
            tau[i, j] = val
            tau[j, i] = val
        if i % 20 == 0:
            print(f"    row {i}/{dim}  ({time.time()-t0:.0f}s)", flush=True)

    print(f"  Built in {time.time()-t0:.1f}s")
    return tau, c_L


def load_zeros(count=200):
    path = os.path.join(DATA_DIR, "zeros_odlyzko_2M")
    zeros = []
    try:
        with open(path) as f:
            for i, line in enumerate(f):
                if i >= count: break
                line = line.strip()
                if line:
                    zeros.append(float(line))
    except:
        path = os.path.join(DATA_DIR, "zeros_100k.txt")
        with open(path) as f:
            for i, line in enumerate(f):
                if i >= count: break
                line = line.strip()
                if line:
                    zeros.append(float(line))
    return np.array(zeros)


def compute_xi_hat(xi_vec, N, L, z):
    """
    ξ̂(z) = ∫_{-L/2}^{L/2} ξ(t) e^{-izt} dt

    In discrete Fourier basis e_j(t) = e^{2πijt/L}/√L:
    ξ̂(z) = Σ_j ξ_j × 2sin(π(j - s))/(2π(j - s)/L)   where s = zL/(2π)

    Using sin(π(j-s)) = (-1)^{j+1} sin(πs) for integer j:
    ξ̂(z) = (√L/π) × sin(πs) × Σ_j (-1)^{j+1} ξ_j/(j-s)   [with alternating signs]

    OR equivalently (standard form from v3):
    ξ̂(z) = (2/√L) × sin(zL/2) × Σ_j ξ_j/(z - 2πj/L)
    """
    pi   = np.pi
    s    = z * L / (2 * pi)
    js   = np.arange(-N, N+1, dtype=float)
    poles = 2*pi*js/L

    # Denominator singularities
    diffs = z - poles  # z - 2πj/L
    mask  = np.abs(diffs) > 1e-10

    # Compute sin(zL/2) prefactor
    sin_prefix = 2.0 / np.sqrt(L) * np.sin(z * L / 2.0)

    if np.abs(np.sin(z*L/2)) < 1e-12:
        # z ≈ 2πk/L (spurious zero) — use L'Hopital
        return 0.0

    val = sin_prefix * np.sum(xi_vec[mask] / diffs[mask])
    return val


def main():
    print("=" * 70)
    print("CONNES Xi COMPARISON: Does ξ̂(z) ≈ Ξ(z)?")
    print("=" * 70)
    print()

    ref_zeros = load_zeros(100)
    print(f"Loaded {len(ref_zeros)} reference zeros")
    print()

    # Test multiple λ values
    configs = [
        {"lambda_sq": 14, "N": 40, "primes": [2,3,5,7,11,13], "label": "λ²=14, N=40"},
        {"lambda_sq": 14, "N": 80, "primes": [2,3,5,7,11,13], "label": "λ²=14, N=80"},
    ]

    all_results = {}

    for cfg in configs:
        lam_sq = cfg["lambda_sq"]
        N      = cfg["N"]
        primes = cfg["primes"]
        label  = cfg["label"]

        L = float(2 * mpmath.log(mpmath.sqrt(lam_sq)))
        print(f"\n{'='*60}")
        print(f"CONFIG: {label}  (L={L:.4f})")
        print(f"{'='*60}")

        tau, c_L = build_weil_matrix(N, L, primes)

        evals, evecs = np.linalg.eigh(tau)
        eps = float(evals[0])
        xi  = evecs[:, 0].copy()

        # Normalize: ⟨1|ξ⟩ = Σ_j ξ_j / √L = 1
        # (Standard Connes normalization: ⟨η|ξ⟩ where η = Σ_j V_j)
        eta_xi = np.sum(xi) / np.sqrt(L)
        if abs(eta_xi) > 1e-12:
            xi /= eta_xi

        print(f"\nε_N = {eps:.8f}")
        print(f"Normalization: ⟨η|ξ⟩ = {eta_xi:.6f} (before norm)")
        print(f"Min eigenvec first 5 components: {xi[:5]}")
        print()

        # ── RATIO TEST: ξ̂(z) / Ξ(z) should be constant ──────────────────
        print("RATIO TEST: ξ̂(z) / Ξ(z)")
        print("(Should be approximately constant if Connes' claim holds)")
        print(f"{'z':>8}  {'ξ̂(z)':>12}  {'Ξ(z)':>12}  {'ratio':>12}")
        print("-" * 52)

        ratios = []
        test_pts = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0]
        for z in test_pts:
            xi_hat_z = compute_xi_hat(xi, N, L, z)
            Xi_z     = float(mpmath.re(Xi(z)))

            if abs(Xi_z) > 1e-20 and abs(xi_hat_z) > 1e-20:
                ratio = xi_hat_z / Xi_z
                ratios.append(ratio)
                print(f"{z:>8.1f}  {xi_hat_z:>12.4e}  {Xi_z:>12.4e}  {ratio:>12.6f}")
            else:
                print(f"{z:>8.1f}  {xi_hat_z:>12.4e}  {Xi_z:>12.4e}  {'N/A':>12}")

        if ratios:
            mean_r = np.mean(ratios)
            std_r  = np.std(ratios)
            cv     = abs(std_r/mean_r) if abs(mean_r) > 1e-12 else float('inf')
            print(f"\nRatio mean={mean_r:.6f}, std={std_r:.6f}, CV={cv:.6f}")
            if cv < 0.01:
                print("*** STRONG: ratio is approximately constant — ξ̂(z) ∝ Ξ(z) ***")
            elif cv < 0.1:
                print("MODERATE: some proportionality detected")
            else:
                print("WEAK: ratio varies significantly — not a good approximation")
        print()

        # ── ZERO COMPARISON ────────────────────────────────────────────────
        print("ZERO COMPARISON: zeros of ξ̂(z) vs zeros of Ξ(z)")
        print("(Scanning z ∈ [12, 200])")

        z_grid = np.linspace(12.0, 200.0, 50000)
        vals = np.array([compute_xi_hat(xi, N, L, z) for z in z_grid])

        found_zeros = []
        for i in range(len(z_grid)-1):
            if vals[i]*vals[i+1] < 0 and not np.isnan(vals[i]) and not np.isnan(vals[i+1]):
                try:
                    z0 = brentq(lambda z: compute_xi_hat(xi, N, L, z),
                               z_grid[i], z_grid[i+1], xtol=1e-10)
                    found_zeros.append(z0)
                except:
                    pass

        print(f"Found {len(found_zeros)} zeros of ξ̂ in [12, 200]")
        print(f"Reference has {len([z for z in ref_zeros if 12<z<200])} zeros in this range")
        print()

        n_tight = 0
        n_loose = 0
        print(f"{'#':>3}  {'Found γ':>10}  {'Ref ζ zero':>10}  {'Error':>10}  Status")
        print("-" * 55)
        for i, (pred, ref) in enumerate(zip(found_zeros[:50], ref_zeros[:50])):
            err = abs(pred - ref)
            if err < 0.1:
                status = "✓✓ (<0.1)"
                n_tight += 1
            elif err < 1.0:
                status = "✓  (<1.0)"
                n_tight += 1
            elif err < 5.0:
                status = "~  (<5.0)"
                n_loose += 1
            else:
                status = "✗  (>5.0)"
            if i < 30:
                print(f"{i+1:>3}  {pred:>10.4f}  {ref:>10.4f}  {err:>10.6f}  {status}")

        print()
        print(f"Tight (err<1.0): {n_tight}/{min(len(found_zeros),50)}")
        print(f"Loose (err<5.0): {n_tight+n_loose}/{min(len(found_zeros),50)}")

        all_results[label] = {
            "N": N, "L": L, "lambda_sq": lam_sq,
            "epsilon_N": eps,
            "n_tight": n_tight,
            "n_loose": n_tight + n_loose,
            "n_found_zeros": len(found_zeros),
            "found_zeros": found_zeros[:60],
            "ratios": ratios,
            "ratio_mean": float(np.mean(ratios)) if ratios else 0,
            "ratio_cv": float(abs(np.std(ratios)/np.mean(ratios))) if ratios and abs(np.mean(ratios)) > 1e-12 else 0,
        }

    # Save
    out = os.path.join(RESULTS_DIR, "connes_xi_compare.json")
    with open(out, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to {out}")

    # ── MPMATH high-precision test ─────────────────────────────────────────
    print()
    print("=" * 70)
    print("HIGH-PRECISION TEST: Can we achieve 10^{-10}+ accuracy on first zero?")
    print("=" * 70)
    print()
    print("The paper claims 10^{-55} error with primes≤13. Testing numerically...")
    print()

    # For λ²=14, N=80: can we refine the first zero to high precision?
    cfg = configs[-1]
    lam_sq = cfg["lambda_sq"]
    N      = cfg["N"]
    L      = float(2 * mpmath.log(mpmath.sqrt(lam_sq)))

    # Re-use last matrix if same config
    print(f"First zeta zero: γ₁ = {ref_zeros[0]:.15f}")
    print(f"First PSWF-based prediction from scan:")

    # Use mpmath to compute Ξ(z) = 0 to high precision (sanity check)
    # Real zero: z = 14.134725141734693...
    z1_ref = float(ref_zeros[0])
    z1_xi  = float(mpmath.re(Xi(z1_ref)))
    print(f"  Ξ(γ₁) = {z1_xi:.6e}  (should be ≈ 0 if γ₁ is correct)")

    # Test if our ξ̂ at λ²=14 vanishes near γ₁
    z_around_first = np.linspace(z1_ref - 2.0, z1_ref + 2.0, 10000)
    # Use result from last config
    xi_last = xi  # last eigenvector computed
    vals_near = [compute_xi_hat(xi_last, N, L, z) for z in z_around_first]
    vals_near = np.array(vals_near)

    has_zero_near = False
    closest_pred = None
    min_err = float('inf')
    for i in range(len(z_around_first)-1):
        if vals_near[i]*vals_near[i+1] < 0:
            has_zero_near = True
            try:
                z0 = brentq(lambda z: compute_xi_hat(xi_last, N, L, z),
                           z_around_first[i], z_around_first[i+1], xtol=1e-13)
                err = abs(z0 - z1_ref)
                if err < min_err:
                    min_err = err
                    closest_pred = z0
            except:
                pass

    if closest_pred is not None:
        print(f"  ξ̂ has a zero near γ₁: {closest_pred:.10f}")
        print(f"  Error: {min_err:.6e}")
        print(f"  Required for paper claim: ~1e-55")
        print(f"  Our precision: ~{min_err:.0e}")
        if min_err < 1.0:
            print("  ✓ Within 1 unit")
        if min_err < 0.01:
            print("  ✓✓ Within 0.01")
    else:
        print("  ✗ No zero of ξ̂ found near γ₁ in [γ₁-2, γ₁+2]")
        print("  This suggests our Weil matrix is not accurately capturing the structure")

    print()
    print("=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print("""
If ratio ξ̂(z)/Ξ(z) ≈ constant:
  → ξ̂ is proportional to Ξ for this λ
  → Strong evidence for Connes' convergence claim
  → Need to prove the ratio → constant as λ→∞ (= Connes' obstacle)

If ratio varies significantly:
  → Our implementation may be incorrect
  → OR the convergence only holds for very large λ
  → OR additional terms are needed

The 10^{-55} accuracy in the paper suggests they use exact analytic methods,
not discrete N-truncation with double precision.
""")

if __name__ == "__main__":
    main()
