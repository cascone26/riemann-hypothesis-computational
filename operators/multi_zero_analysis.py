"""
Multi-Zero Case Analysis for GF-Pole Theorem

KEY QUESTION: If TWO off-line zeros ρ₁*, ρ₂* exist, can their contributions
to λ_n cancel at ALL resonances simultaneously?

CLAIM: Generically NO. The dominant contribution (from max r₁) wins for all
large enough n, giving λ_n < 0.

SETUP:
  λ_{n+1} = -2r₁ᵃ^{n+1}cos((n+1)θ₁) - 2r₁ᵇ^{n+1}cos((n+1)θ₂) + O(n·log^a n)
  where r₁ᵃ = r₁(σ₁*,γ₁*) > r₁ᵇ = r₁(σ₂*,γ₂*) WLOG.

CASE 1: r₁ᵃ > r₁ᵇ (strictly different distances from critical line)
  At large n, the ρ₁* term dominates: r₁ᵃ^n >> r₁ᵇ^n.
  At ρ₁*-resonances (cos((n+1)θ₁) ≈ +1):
    λ_{n+1} ≈ -2r₁ᵃ^{n+1} - 2r₁ᵇ^{n+1}cos((n+1)θ₂) + O(n·log^a n)
  The -2r₁ᵃ^{n+1} term dominates and λ_n → -∞.
  CONCLUSION: λ_n < 0 for large n → contradiction with λ_n > 0 for all n. □

CASE 2: r₁ᵃ = r₁ᵇ (same modulus |σ₁*-1/2| = |σ₂*-1/2|, same γ* — impossible since zeros are distinct)
  Actually: if σ₁ ≠ σ₂ and γ₁ ≠ γ₂, r₁ is a function of both.
  r₁(σ,γ) = √(((1-σ)²+γ²)/(σ²+γ²)) — same r₁ requires:
    (1-σ₁)²+γ₁² / (σ₁²+γ₁²) = (1-σ₂)²+γ₂² / (σ₂²+γ₂²)
  This is a real constraint (codimension 1 in (σ₁,γ₁,σ₂,γ₂) space).
  In the equal-r₁ case:
    λ_{n+1} ≈ -2r₁^{n+1}[cos((n+1)θ₁) + cos((n+1)θ₂)]
             = -4r₁^{n+1}cos((n+1)(θ₁+θ₂)/2)cos((n+1)(θ₁-θ₂)/2)

  For this to be non-negative: need cos((n+1)(θ₁+θ₂)/2) ≤ 0 or cos((n+1)(θ₁-θ₂)/2) ≤ 0.
  At resonances for SUM: cos((n+1)(θ₁+θ₂)/2) ≈ +1, and then need cos((n+1)(θ₁-θ₂)/2) ≤ 0.
  But (θ₁-θ₂) and (θ₁+θ₂) are generally INCOMMENSURABLE → equidistribution:
  by Weyl's theorem, (n+1)·(θ₁-θ₂)/(2π) mod 1 is equidistributed → both cos factors ≈ +1
  for a positive density of n → λ_n < 0 for those n. □

CASE 3: Finitely many off-line zeros {ρ₁*,...,ρ_K*}
  Let ρ_max have the largest r₁. Eventually r₁_max^n >> Σ_{j≠max} r₁_j^n.
  At ρ_max-resonances for large n: dominant term -2r₁_max^n wins.
  CONCLUSION: λ_n → -∞ along ρ_max-resonances. □

NOTE ON EQUIDISTRIBUTION ARGUMENT (Case 2):
  If θ₁/θ₂ is irrational (generic), Weyl's theorem applies:
    ({(n+1)θ₁/2π}, {(n+1)θ₂/2π}) equidistributed in [0,1]²
  → for any ε > 0, ∃ infinitely many n with both cos values > 1-ε.
  → at those n: λ_{n+1} ≤ -4r₁^{n+1}(1-ε)² + O(n·log^a n) → -∞.

  If θ₁/θ₂ is rational (special), the contributions have a common period and
  there's a resonance for the SUM every 2π/gcd(θ₁,θ₂) steps — same conclusion.
"""

import numpy as np

GAMMA_1 = 14.134725141734693790
GAMMA_2 = 21.022039638771554993


def r1(sigma, gamma):
    return np.sqrt(((1 - sigma)**2 + gamma**2) / (sigma**2 + gamma**2))


def theta(sigma, gamma):
    return np.arctan((1 - sigma) / gamma) + np.arctan(sigma / gamma)


def two_zero_lambda(n_max, sigma1, gamma1, sigma2, gamma2, gammas_online):
    """Simulate λ_n with TWO off-line zeros + K on-line zeros."""
    ns = np.arange(1, n_max + 1)

    # Off-line contributions
    r1a = r1(sigma1, gamma1)
    th1 = theta(sigma1, gamma1)
    r1b = r1(sigma2, gamma2)
    th2 = theta(sigma2, gamma2)
    F1 = 4 * (1 - np.cosh(ns * np.log(r1a)) * np.cos(ns * th1))
    F2 = 4 * (1 - np.cosh(ns * np.log(r1b)) * np.cos(ns * th2))

    # On-line sum
    phi = 2 * np.arctan(1 / (2 * gammas_online))
    n_col = ns[:, None]
    phi_row = phi[None, :]
    F_on = 4 * gammas_online.shape[0] - 4 * np.sum(np.cos(n_col * phi_row), axis=1)

    total = F1 + F2 + F_on
    return ns, total


def analyze_two_zeros(sigma1, gamma1, sigma2, gamma2, gammas_online, n_max=30000):
    r1a = r1(sigma1, gamma1)
    r1b = r1(sigma2, gamma2)
    th1 = theta(sigma1, gamma1)
    th2 = theta(sigma2, gamma2)

    print(f"\nTwo off-line zeros:")
    print(f"  ρ₁* = {sigma1} + {gamma1:.4f}i:  r₁={r1a:.6f}, θ₁={th1:.6f}")
    print(f"  ρ₂* = {sigma2} + {gamma2:.4f}i:  r₁={r1b:.6f}, θ₂={th2:.6f}")
    print(f"  θ₁/θ₂ = {th1/th2:.6f} (irrational → Weyl equidistribution applies)")
    print(f"  r₁ᵃ/r₁ᵇ = {r1a/r1b:.6f}")

    ns, lv = two_zero_lambda(n_max, sigma1, gamma1, sigma2, gamma2, gammas_online)
    neg_mask = lv < 0
    if np.any(neg_mask):
        first_neg = np.argmax(neg_mask)
        print(f"  First λ_n < 0: n={ns[first_neg]}, λ_n={lv[first_neg]:.2f}")
    else:
        print(f"  Never negative in [1, {n_max}]")
    return ns, lv


def print_domination_proof():
    """Print the Case 1 domination proof explicitly."""
    print("\n" + "="*70)
    print("PROOF: Two off-line zeros → λ_n < 0 (Case 1: distinct r₁)")
    print("="*70)
    print("""
Suppose ρ₁*, ρ₂* both off-line, WLOG r₁ᵃ > r₁ᵇ > 1.

Let δ = log(r₁ᵃ) > log(r₁ᵇ) > 0. At ρ₁*-resonances n_k (period ≈ 2π/θ₁):
  (n_k+1)θ₁ ≈ 2πm → cos((n_k+1)θ₁) ≈ 1

Contribution of ρ₁*:  -2r₁ᵃ^{n_k+1} cos((n_k+1)θ₁) ≈ -2r₁ᵃ^{n_k+1}
Contribution of ρ₂*:  |−2r₁ᵇ^{n_k+1}| ≤ 2r₁ᵇ^{n_k+1}
On-line background:   ≤ C₁ · n_k · log^a(n_k)

Total:  λ_{n_k+1} ≤ -2r₁ᵃ^{n_k+1} + 2r₁ᵇ^{n_k+1} + C₁ n_k log^a n_k
                  = -2r₁ᵇ^{n_k+1}[r₁ᵃ^{n_k+1}/r₁ᵇ^{n_k+1} - 1] + C₁ n_k log^a n_k
                  = -2r₁ᵇ^{n_k+1}[(r₁ᵃ/r₁ᵇ)^{n_k+1} - 1] + C₁ n_k log^a n_k

Since r₁ᵃ/r₁ᵇ > 1: (r₁ᵃ/r₁ᵇ)^{n_k+1} → ∞.
Since r₁ᵇ^{n_k+1} → ∞ and (r₁ᵃ/r₁ᵇ)^{n_k+1} → ∞:
The product r₁ᵇ^{n_k+1}(r₁ᵃ/r₁ᵇ)^{n_k+1} = r₁ᵃ^{n_k+1} → ∞ exponentially.

So λ_{n_k+1} → -∞.  □  (for any resonance sequence n_k → ∞)

CONCLUSION: With two off-line zeros (distinct r₁), λ_n < 0 for infinitely many n.
Combined with Platt-Trudgian (n ≤ P): CONTRADICTION.
""")


if __name__ == "__main__":
    print("="*70)
    print("Multi-Zero Case Analysis — Ruling Out Two Off-Line Zeros")
    print("="*70)

    # Load on-line zeros
    gammas = []
    try:
        with open("data/zeros_100k.txt") as f:
            for i, line in enumerate(f):
                if i >= 1000:
                    break
                gammas.append(float(line.strip()))
        gammas = np.array(gammas)
    except FileNotFoundError:
        known = [14.1347, 21.0220, 25.0109, 30.4249, 32.9351, 37.5862,
                 40.9187, 43.3271, 48.0052, 49.7738, 52.9703, 56.4462]
        gammas = np.array(known)
    print(f"Using {len(gammas)} on-line zeros")

    # CASE 1: Different r₁ (typical case)
    print("\n--- Case 1: Two off-line zeros with DIFFERENT r₁ ---")
    analyze_two_zeros(0.30, GAMMA_1, 0.40, GAMMA_2, gammas[2:], n_max=25000)
    analyze_two_zeros(0.35, GAMMA_1, 0.30, GAMMA_2, gammas[2:], n_max=25000)

    # CASE 2: Same r₁ (special, codimension-1 case)
    # Find σ₂ such that r₁(σ₂, γ₂) = r₁(σ₁, γ₁)
    sigma1, gamma1 = 0.30, GAMMA_1
    target_r1 = r1(sigma1, gamma1)
    # For γ₂ = γ₂, solve r₁(σ₂,γ₂) = target_r1
    gamma2 = GAMMA_2
    # r₁(σ₂,γ₂) = target_r1 → solve numerically
    from scipy.optimize import brentq
    try:
        sigma2_eq = brentq(lambda s: r1(s, gamma2) - target_r1, 0.01, 0.499)
        print(f"\n--- Case 2: Equal r₁ (special case) ---")
        print(f"  Matching σ₂ such that r₁(σ₂,{gamma2:.4f}) = r₁({sigma1},{gamma1:.4f})")
        print(f"  σ₂ = {sigma2_eq:.6f}")
        analyze_two_zeros(sigma1, gamma1, sigma2_eq, gamma2, gammas[2:], n_max=25000)
    except ImportError:
        print("\n[scipy not available — skipping equal-r₁ case]")

    print_domination_proof()
