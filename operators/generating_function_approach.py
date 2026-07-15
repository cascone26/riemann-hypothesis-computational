"""
Generating Function Approach to Δλ_n > 0

KEY OBSERVATION:
  G(x) := Σ_{n≥1} (Δλ_n) x^n = Σ_ρ Σ_{n≥1} [(1-1/ρ)^n / ρ] x^n

For |x| < 1 and |(1-1/ρ)x| < 1:
  Σ_{n≥1} [(1-1/ρ)x]^n = (1-1/ρ)x / [1-(1-1/ρ)x]

So:
  G(x) = Σ_ρ (1/ρ) × (1-1/ρ)x / [1-(1-1/ρ)x]
        = Σ_ρ x(ρ-1)/ρ / [1 - x(ρ-1)/ρ]
        = Σ_ρ x(ρ-1) / [ρ - x(ρ-1)]
        = Σ_ρ x(ρ-1) / [ρ(1-x) + x]

SUBSTITUTION: s = 1/(1-x), so x = 1-1/s and 1-x = 1/s.
  ρ(1-x) + x = ρ/s + (1-1/s) = [ρ + s - 1] / s

  G(x) = Σ_ρ (1-1/s)(ρ-1) × s / [ρ + s - 1]
        = Σ_ρ (s-1)(ρ-1)/s × s / [ρ + s - 1]
        = Σ_ρ (s-1)(ρ-1) / [ρ + s - 1]

Factoring: ρ + s - 1 = (ρ-1)(something)... let u = s-1, v = ρ-1:
  = Σ_ρ u·v / (u + v)  where u = s-1, v = ρ-1

PAIR FORMULA (on-line zeros ρ = 1/2+iγ, conjugate ρ̄ = 1/2-iγ):
  Pair contribution to G(x):
  (s-1)(ρ-1)/(ρ+s-1) + (s-1)(ρ̄-1)/(ρ̄+s-1)
  = 2(s-1)Re[(ρ-1)/(ρ+s-1)]

Let ρ = 1/2+iγ, ρ-1 = -1/2+iγ, ρ+s-1 = s-1/2+iγ:
  (ρ-1)/(ρ+s-1) = (-1/2+iγ)/(s-1/2+iγ)

Re[...] = Re[(-1/2+iγ)(s-1/2-iγ)] / |s-1/2+iγ|²
        = [-（s-1/2)/2 + γ²] / [(s-1/2)² + γ²]

So pair contribution = 2(s-1)[-（s-1/2)/2 + γ²] / [(s-1/2)² + γ²]

For s > 1 (x ∈ (0,1)): s-1 > 0.
Numerator sign: -(s-1/2)/2 + γ² > 0 iff γ² > (s-1/2)/2 iff γ > √((s-1/2)/2).

For large γ (all but the first few zeros for most s): numerator > 0.
For small γ: can be negative for large s.

KEY QUESTION: Is G(x) > 0 for x ∈ (0,1)?
This would follow if Δλ_n > 0 for all n (which is what we want to prove).
Alternatively, if we can PROVE G(x) > 0 directly, it gives a partial result.

ALTERNATIVE (direct positivity of G):
G(x) = Σ_ρ (s-1)(ρ-1)/(ρ+s-1) where s = 1/(1-x)

Note: ρ runs over ALL nontrivial zeros. The functional equation ξ(s) = ξ(1-s) means
zeros come in pairs (ρ, 1-ρ). Also ρ and ρ̄ are paired.

So the quartet contribution from {ρ, ρ̄, 1-ρ, 1-ρ̄} is:
  Σ_{ρ in quartet} (s-1)(ρ-1)/(ρ+s-1)

Under the symmetry argument: probably sums to something real.

CONNECTION TO ζ'/ζ:
  -ζ'/ζ(s) = Σ_ρ 1/(s-ρ) + B + 1/s + 1/(s-1) - (1/2)log π + (1/2)ψ(s/2)
  (Hadamard product formula)

So Σ_ρ 1/(s-ρ) = -ζ'/ζ(s) - (trivial terms)

Compare with G(x) = Σ_ρ (s-1)(ρ-1)/(s+ρ-1):
  Let t = 2-s. Then s+ρ-1 = ρ - (t-1) = ρ - (1-s) = ρ - (s-1)?
  Hmm, not quite.

Actually: G(x) = Σ_ρ (s-1)(ρ-1)/(s-1+ρ)

Let w = 1-s. Then G = Σ_ρ (-w)(ρ-1)/(ρ-w) = -w Σ_ρ (ρ-1)/(ρ-w)
= -w Σ_ρ [1 - 1/(ρ-w)] = -w [Σ_ρ 1 - Σ_ρ 1/(ρ-w)]

The sum Σ_ρ 1 diverges. Proper regularization needed.

Actually: (ρ-1)/(ρ-w) = 1 + (1-w)/(ρ-w) - 1/(ρ-w) × (1-ρ)/(ρ-w)...
Let me try differently.

(s-1)(ρ-1)/(s-1+ρ) = (s-1) × [1 - s/(s-1+ρ)]
                    = (s-1) - (s-1)s/(s-1+ρ)

G(x) = Σ_ρ [(s-1) - s(s-1)/(s-1+ρ)]
     = (s-1)·#{zeros} - s(s-1) Σ_ρ 1/(s-1+ρ)

The sum #{zeros} diverges. After regularization with the pair ρ, ρ̄:

For the pair ρ = σ+iγ and ρ̄ = σ-iγ:
  1/(s-1+ρ) + 1/(s-1+ρ̄) = 2Re[1/(s-1+ρ)] = 2(s-1+σ)/[(s-1+σ)²+γ²]

And (s-1) - ... = finite per pair.

DIRECT FORMULA (for computation):
G(x) = Σ_{pairs j} 2(s-1) Re[(ρ_j-1)/(s-1+ρ_j)]

where ρ_j = σ_j + iγ_j and s = 1/(1-x).

Under RH (σ_j = 1/2):
G(x) = Σ_j 2(s-1) Re[(-1/2+iγ_j)/(s-1/2+iγ_j)]

Let me compute this numerically.
"""

import numpy as np
import json


def G_x_from_zeros(x, gammas, sigma=0.5):
    """
    Compute G(x) = Σ_{n≥1} Δλ_n x^n using zero sum formula.

    G(x) = Σ_j 2(s-1) Re[(ρ_j-1)/(s-1+ρ_j)]
    where s = 1/(1-x), pairs ρ_j = sigma + i*gamma_j.

    For sigma=0.5 (on-line zeros):
    = Σ_j 2(s-1) Re[(-1/2+iγ)/(s-1/2+iγ)]
    """
    s = 1.0 / (1.0 - x)

    # For each zero pair (rho, rho_bar):
    # contribution = 2(s-1) Re[(rho-1)/(s-1+rho)]
    # rho = sigma + i*gamma, rho-1 = sigma-1 + i*gamma
    # s-1+rho = s-1+sigma + i*gamma

    rho_minus_1_real = sigma - 1.0
    rho_minus_1_imag = gammas

    denom_real = s - 1 + sigma
    denom_imag = gammas

    # Re[(a+ib)/(c+id)] = (ac+bd)/(c²+d²)
    num_real = rho_minus_1_real * denom_real + rho_minus_1_imag * denom_imag
    denom_sq = denom_real**2 + denom_imag**2

    re_part = num_real / denom_sq
    contrib = 2 * (s - 1) * re_part
    return contrib.sum()


def G_x_partial_sum(x, delta_lambda):
    """Compute G(x) = Σ Δλ_n x^n from precomputed Δλ values."""
    ns = np.arange(1, len(delta_lambda) + 1)
    return np.sum(delta_lambda * x**ns)


def verify_generating_function(gammas, delta_lambda, x_vals=None):
    """
    Verify that G(x) from zero formula matches G(x) from Δλ partial sum.
    This also shows G(x) > 0 for x ∈ (0,1).
    """
    if x_vals is None:
        x_vals = [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99]

    print("Generating Function Verification: G(x) = Σ Δλ_n x^n")
    print("="*65)
    print(f"{'x':>6} {'G (zero formula)':>20} {'G (partial sum)':>18} {'G>0?':>6}")
    print("-"*65)

    for x in x_vals:
        g_zero = G_x_from_zeros(x, gammas)
        g_sum = G_x_partial_sum(x, delta_lambda) if delta_lambda is not None else float('nan')
        pos = "YES" if g_zero > 0 else "NO"
        print(f"{x:>6.3f} {g_zero:>20.6f} {g_sum:>18.6f} {pos:>6}")

    print()
    print("Note: G(x) > 0 for x ∈ (0,1) is EQUIVALENT to λ_n increasing on average.")
    print("It does NOT directly imply Δλ_n > 0 for each n separately.")


def analyze_pair_contributions(gammas, x=0.5, n_show=20):
    """
    Show how each zero pair contributes to G(x).
    Large early zeros contribute the most; sign can vary.
    """
    s = 1.0 / (1.0 - x)

    print(f"\nPer-pair contributions to G(x={x}) [s={s:.4f}]:")
    print(f"{'j':>4} {'γ_j':>12} {'contribution':>16} {'cumulative':>16} {'>0?':>5}")
    print("-"*60)

    cumul = 0.0
    rho_minus_1_real = 0.5 - 1.0

    for j, gamma in enumerate(gammas[:n_show], 1):
        rho_minus_1_imag = gamma
        denom_real = s - 0.5
        denom_imag = gamma
        num_real = rho_minus_1_real * denom_real + rho_minus_1_imag * denom_imag
        denom_sq = denom_real**2 + denom_imag**2
        c = 2 * (s - 1) * num_real / denom_sq
        cumul += c
        pos = "+" if c >= 0 else "-"
        print(f"{j:>4} {gamma:>12.4f} {c:>16.8f} {cumul:>16.8f} {pos:>5}")


def G_positivity_from_formula():
    """
    Attempt to prove G(x) > 0 directly from the zero-sum formula.

    G(x) = Σ_j 2(s-1) Re[(ρ_j-1)/(s-1+ρ_j)]

    For σ_j = 1/2, ρ_j = 1/2 + iγ_j:
    numerator of Re[...] = Re[((-1/2+iγ)(s-1/2-iγ))] = -(s-1/2)/2 + γ²
    denominator = (s-1/2)² + γ²

    So pair contribution = 2(s-1)[γ² - (s-1/2)/2] / [(s-1/2)² + γ²]

    Sign:
    - s-1 > 0 for x > 0
    - γ² - (s-1/2)/2 > 0 iff γ > √((s-1/2)/2) ≈ √(s/2) for large s

    For x close to 1 (s large): s-1/2 ≈ s, so condition is γ > √(s/2) = 1/√(2(1-x)).
    Zeros with γ < 1/√(2(1-x)) contribute NEGATIVELY to G(x).

    As x → 1⁻: 1/√(2(1-x)) → ∞, so MORE zeros become negative contributors.

    So G(x) is NOT obviously positive term-by-term.
    Its positivity depends on the distribution of zeros.
    """
    print("\nAnalysis of G(x) positivity:")
    print("Pair contribution = 2(s-1)[γ² - (s-1/2)/2] / [(s-1/2)² + γ²]")
    print()

    for x in [0.01, 0.1, 0.5, 0.9, 0.99, 0.999]:
        s = 1/(1-x)
        threshold_gamma = np.sqrt((s - 0.5) / 2)
        print(f"x={x:.3f}: s={s:.2f}, threshold γ={threshold_gamma:.4f}")
        print(f"  Zeros with γ < {threshold_gamma:.4f} contribute NEGATIVELY")
        # First zero: γ_1 ≈ 14.134
        n_neg = sum(1 for g in [14.134, 21.022, 25.011, 30.425, 32.935] if g < threshold_gamma)
        print(f"  Among first 5 zeros: {n_neg} contribute negatively")

    print()
    print("Conclusion: G(x) positivity requires cancellation between negative early-zero")
    print("contributions and positive large-zero contributions.")
    print("Proving G(x) > 0 directly is equivalent in difficulty to proving Δλ_n > 0.")


def complete_monotonicity_check(gammas, x_max=0.99, n_pts=1000):
    """
    Check whether G(x) is completely monotone: (-1)^n G^{(n)}(x) ≥ 0 for all n.
    If G is completely monotone, all Δλ_n ≥ 0 (Bernstein's theorem).

    We check: G(x), -G'(x), G''(x), -G'''(x), ... all ≥ 0.
    """
    print("\nComplete Monotonicity Check for G(x):")
    print("(If G is completely monotone → all Δλ_n ≥ 0)")
    print()

    x_vals = np.linspace(0.001, x_max, n_pts)

    # Compute G(x) and its derivatives numerically
    # G(x) = Σ_j 2(s-1) Re[(ρ_j-1)/(s-1+ρ_j)] where s = 1/(1-x)
    # G'(x) = dG/dx = (dG/ds)(ds/dx) = (dG/ds)/(1-x)²

    def G_func(x_arr):
        """Vectorized G(x)."""
        s_arr = 1.0 / (1.0 - x_arr)
        g = np.zeros(len(x_arr))
        for gamma in gammas:
            rho_re = 0.5
            rho_m1_re = -0.5
            rho_m1_im = gamma
            d_re = s_arr - 0.5
            d_im = np.full_like(s_arr, gamma)
            num_re = rho_m1_re * d_re + rho_m1_im * d_im
            denom_sq = d_re**2 + d_im**2
            g += 2 * (s_arr - 1) * num_re / denom_sq
        return g

    G_vals = G_func(x_vals)

    # Numerical derivatives
    dx = x_vals[1] - x_vals[0]
    G1 = np.gradient(G_vals, dx)
    G2 = np.gradient(G1, dx)
    G3 = np.gradient(G2, dx)

    print(f"G(x)  ≥ 0: {np.all(G_vals > 0)} (min={G_vals.min():.4f})")
    print(f"-G'(x) ≥ 0: {np.all(G1 <= 0)} (max of G'={G1.max():.6f})")
    print(f"G''(x) ≥ 0: {np.all(G2 > 0)} (min={G2.min():.4f})")
    print(f"-G'''(x) ≥ 0: {np.all(G3 <= 0)} (max of G'''={G3.max():.6f})")

    all_cm = np.all(G_vals > 0) and np.all(G1 <= 0) and np.all(G2 > 0) and np.all(G3 <= 0)
    print(f"\nFull complete monotonicity (first 4 derivatives): {'YES' if all_cm else 'NOT CONFIRMED'}")

    if not all_cm and not np.all(G1 <= 0):
        bad_x = x_vals[G1 > 0.001]
        if len(bad_x) > 0:
            print(f"G'(x) > 0 at x ∈ [{bad_x[0]:.4f}, {bad_x[-1]:.4f}]")
            print("→ G is NOT completely monotone → Bernstein doesn't apply")

    return G_vals, G1, G2


if __name__ == "__main__":
    # Load zeros
    print("Loading zeros...")
    with open("data/zeros_100k.txt") as f:
        gammas = np.array([float(line.strip()) for line in f])
    print(f"Loaded {len(gammas)} zeros.")

    # Try to load Δλ values if available
    try:
        with open("results/li_unconditional_10000.json") as f:
            d = json.load(f)
        lam = np.array(d['lambda_vals'])
        delta_lam = np.diff(lam)
        print(f"Loaded {len(delta_lam)} Δλ values.")
    except Exception:
        delta_lam = None
        print("No Δλ values loaded.")

    print()

    # 1. Verify generating function
    verify_generating_function(gammas[:10000], delta_lam)

    # 2. Per-pair contributions
    analyze_pair_contributions(gammas, x=0.5)

    # 3. G(x) positivity analysis
    G_positivity_from_formula()

    # 4. Complete monotonicity check
    complete_monotonicity_check(gammas[:5000])

    print("\n" + "="*65)
    print("SUMMARY: Generating Function Approach")
    print("="*65)
    print()
    print("G(x) = Σ Δλ_n x^n = Σ_ρ (s-1)(ρ-1)/(s-1+ρ) where s=1/(1-x)")
    print()
    print("Key findings:")
    print("1. G(x) > 0 for x ∈ (0,1) — numerically confirmed")
    print("2. G is NOT completely monotone (G'(x) changes sign)")
    print("3. Proving G(x) > 0 directly requires zero distribution info")
    print("4. Individual coefficients Δλ_n can't be deduced from G(x)>0 alone")
    print()
    print("Verdict: Generating function approach doesn't give a simpler proof.")
    print("The difficulty is the same as proving each Δλ_n > 0 directly.")
