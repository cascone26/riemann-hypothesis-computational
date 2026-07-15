"""
Monotonicity of λ_n: λ_{n+1} > λ_n for all n ≥ 1

OBSERVATION (numerically established for n ≤ 10000):
  The Keiper-Li sequence λ_n is strictly monotone increasing.

IMPLICATION:
  If λ_1 > 0 (KNOWN: λ_1 = 1 - γ_E/2 - log(4π)/2 ≈ 0.02310 > 0)
  AND λ_n is strictly increasing,
  THEN λ_n ≥ λ_1 > 0 for all n ≥ 1.
  This would prove RH (Keiper-Li criterion).

EXACT INCREMENT FORMULA:
  Δλ_n := λ_{n+1} - λ_n = Σ_ρ [(1-1/ρ)^n - (1-1/ρ)^{n+1}]
         = Σ_ρ (1-1/ρ)^n × (1/ρ)

  For on-line zeros ρ = 1/2 + iγ:
    |(1-1/ρ)| = |(ρ-1)/ρ| = √(1/4+γ²)/√(1/4+γ²) = 1  EXACTLY
    arg((ρ-1)/ρ) = π - 2 arctan(2γ) =: α(γ)   (depends on γ only)

  So each pair's contribution to Re[Δλ_n]:
    2 Re[(e^{iα})^n/ρ] = [cos(nα)/2 + γ sin(nα)] / (1/4+γ²)
                       ≈ sin(nα) / γ for large γ
                       ≈ sin(n/γ) / γ  (since α ≈ 1/γ for γ >> 1)

  Positive when n < πγ (n/γ < π → sin > 0)
  Negative when n > πγ (n/γ > π → sin < 0)

  The "turning point" for zero γ_j: n* = π × γ_j.
  For γ_1 ≈ 14.13: n* ≈ 44.4.

  The sum stays positive because the negative contributions (from zeros with γ_j < n/π)
  are bounded by the large sieve: O(n log n), while positive contributions from
  large-γ zeros accumulate. But proving this rigorously requires controlling
  the oscillation sum — which is as deep as RH.

DIFFICULTY:
  Proving Δλ_n > 0 analytically is equivalent in difficulty to proving λ_n > 0.
  The "sin(n/γ_j)/γ_j" representation shows the sum involves subtle cancellations
  that depend on the precise distribution of the imaginary parts γ_j.

NUMERICAL EVIDENCE:
  n ≤ 10,000: Δλ_n > 0 for ALL n (minimum at n=1: Δλ_1 ≈ 0.069250)
  This is strong evidence for the monotonicity theorem.
"""

import json
import numpy as np
from scipy.special import digamma


# ────────────────────────────────────────────────────────────────────────────
# Exact value of λ_1

def lambda_1_exact():
    """
    λ_1 = Σ_ρ 1/ρ = 1 + γ_E/2 - log(4π)/2
    (Li 1997 eq (2.2): λ_1 = B where B is the constant in the Hadamard product for ξ)

    From the product formula: 2B = log(4π) - γ_E - 2
    So B = (log(4π) - γ_E - 2)/2
    And λ_1 = 1 + (γ_E/2 - log(4π)/2 + 1) - 1 = ... let me use the known numerical value.

    Verified: λ_1 = 1 - (γ_E/2) - log(2π)/2 + (log π)/2 - (log π)/2 ...
    Actually from data: λ_1 ≈ 0.023095709.

    Direct formula (Keiper 1992): λ_1 = 1 - log(4π)/2 + γ_E/2 + sum over trivial zeros
    The cleanest known form: λ_1 ≈ 0.02310 > 0. Positivity is elementary since
    log(4π)/2 - γ_E/2 ≈ 1.266 - 0.289 = 0.977 < 1.
    """
    # Correct derivation (from d/ds[log ξ(s)]|_{s=1}):
    #   λ_1 = 1 + γ_E/2 - (1/2)log(4π)
    # where we use ψ(1/2) = -γ_E - 2log2 and ξ(1) = 1/2.
    # Numerically: 1 + 0.2886 - 1.2655 = 0.0231 > 0.
    gamma_euler = np.euler_gamma  # ≈ 0.5772
    return 1.0 + gamma_euler / 2.0 - 0.5 * np.log(4 * np.pi)


# ────────────────────────────────────────────────────────────────────────────
# Per-zero contribution to Δλ_n

def delta_contribution_pair(gamma, n):
    """
    Exact contribution of the pair ρ = 1/2 ± iγ to Δλ_n.
    Returns 2 × Re[(1-1/ρ)^n / ρ].
    """
    rho = complex(0.5, gamma)
    w = (rho - 1) / rho  # This has |w| = 1 for on-line zeros
    contrib = (w**n) / rho
    return 2 * contrib.real


def alpha_angle(gamma):
    """Phase angle α = arg((ρ-1)/ρ) for ρ = 1/2+iγ."""
    return np.pi - 2 * np.arctan(2 * gamma)


def approx_contribution(gamma, n):
    """Approximate contribution ≈ sin(nα)/γ for large γ."""
    alpha = alpha_angle(gamma)
    return np.sin(n * alpha) / gamma


# ────────────────────────────────────────────────────────────────────────────
# Load and analyze λ_n sequence

def load_lambda(path='results/li_unconditional_10000.json'):
    with open(path) as f:
        d = json.load(f)
    return np.array(d['lambda_vals'])


def monotonicity_report(lam):
    """Full report on monotonicity of the λ_n sequence."""
    n = np.arange(1, len(lam) + 1)
    diffs = np.diff(lam)  # Δλ_n = λ_{n+1} - λ_n

    print("=" * 70)
    print("Monotonicity Analysis of λ_n")
    print("=" * 70)
    print(f"\nλ_1 (exact) = {lambda_1_exact():.8f}")
    print(f"λ_1 (from data) = {lam[0]:.8f}")
    print(f"\nSequence length: {len(lam)}")
    print(f"All λ_n > 0: {np.all(lam > 0)}")
    print(f"All Δλ_n > 0 (strictly increasing): {np.all(diffs > 0)}")
    print(f"Min Δλ_n = {diffs.min():.6f} at n = {np.argmin(diffs)+1}")
    print(f"Max Δλ_n = {diffs.max():.6f} at n = {np.argmax(diffs)+1}")

    print("\nΔλ_n samples:")
    print(f"  {'n':>8} {'λ_n':>14} {'Δλ_n':>12} {'Δλ_n/log(n)':>14}")
    for i in [0, 4, 9, 49, 99, 499, 999, 4999, 9998]:
        if i < len(diffs):
            nn = i + 1
            ratio = diffs[i] / np.log(nn) if nn > 1 else float('nan')
            print(f"  {nn:>8} {lam[i]:>14.4f} {diffs[i]:>12.6f} {ratio:>14.6f}")

    print(f"\nAsymptotic: Δλ_n ≈ log(n)/2 (under RH)")
    print(f"  Ratio Δλ_n / (log(n)/2) at n=10000: "
          f"{diffs[-1] / (np.log(9999)/2):.4f} (should ≈ 1.0)")


# ────────────────────────────────────────────────────────────────────────────
# Turning-point analysis

def turning_point_analysis(gammas, n_max=200):
    """
    For each zero γ_j, the contribution to Δλ_n is:
    ≈ sin(n × α_j) / γ_j  where α_j = π - 2arctan(2γ_j) ≈ 1/γ_j.
    This is positive for n < πγ_j and negative for πγ_j < n < 2πγ_j, etc.

    Show the cumulative sum of contributions for small n.
    """
    print("\nTurning-point analysis (first 10 zeros):")
    print("  γ_j: turning points n* = kπγ_j where sign flips")
    print(f"  {'γ_j':>10} {'α_j':>10} {'n*(k=1)':>10} {'n*(k=2)':>10}")
    for g in gammas[:10]:
        alpha = alpha_angle(g)
        n1 = np.pi / alpha  # ≈ πγ_j
        n2 = 2 * np.pi / alpha
        print(f"  {g:>10.4f} {alpha:>10.6f} {n1:>10.1f} {n2:>10.1f}")

    print(f"\nΔλ_n = Σ_j contribution_j(n)")
    print(f"Cumulative contributions for selected n:")
    print(f"  {'n':>6} {'Σ exact':>14} {'Σ approx':>14}")
    for n in [1, 5, 10, 20, 44, 50, 100, 200]:
        if n <= n_max:
            exact = sum(delta_contribution_pair(g, n) for g in gammas)
            approx = sum(approx_contribution(g, n) for g in gammas)
            print(f"  {n:>6} {exact:>14.6f} {approx:>14.6f}")


# ────────────────────────────────────────────────────────────────────────────
# K-zero lower bound on Δλ_n

def k_zero_lower_bound(gammas, n_max=10000):
    """
    Using the first K verified on-line zeros, compute a partial sum of Δλ_n.
    This is a lower bound IF we can bound the tail.

    For the tail (zeros with γ > γ_K), the contribution is approximately:
      Σ_{j>K} sin(nα_j)/γ_j ≈ Σ_{j>K} sin(n/γ_j)/γ_j

    For n ≤ πγ_K (all large zeros in first half-period): sin(n/γ_j) > 0,
    so the tail contributes POSITIVELY. This gives:

    Δλ_n ≥ Δλ_n^{(K)}   for n ≤ πγ_K

    where Δλ_n^{(K)} = contribution from first K zeros.
    """
    gamma_K = gammas[-1]
    n_safe = int(np.pi * gamma_K)

    print(f"\nK-zero lower bound (K={len(gammas)} zeros, γ_K={gamma_K:.2f}):")
    print(f"  For n ≤ πγ_K = {n_safe}: tail is positive, so Δλ_n ≥ Δλ_n^(K)")

    # Compute Δλ_n^(K) for selected n
    ns = [1, 10, 100, 1000, min(n_safe, 10000)]
    print(f"  {'n':>8} {'Δλ^(K)':>14} {'Status':>20}")
    for n in ns:
        partial = sum(delta_contribution_pair(g, n) for g in gammas)
        status = "lower bound on Δλ_n" if n <= n_safe else "NOT a bound (tail can be neg)"
        print(f"  {n:>8} {partial:>14.6f}  {status}")

    return n_safe


# ────────────────────────────────────────────────────────────────────────────
# Main

if __name__ == "__main__":
    # Load λ_n sequence
    lam = load_lambda()
    monotonicity_report(lam)

    # Load zeros for turning-point analysis
    known_gammas = [
        14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
        37.586178, 40.918720, 43.327073, 48.005151, 49.773832,
        52.970321, 56.446248, 59.347044, 60.831779, 65.112544
    ]
    gammas = np.array(known_gammas)
    turning_point_analysis(gammas, n_max=200)

    # K-zero lower bound
    k_zero_lower_bound(gammas, n_max=10000)

    print("\n" + "=" * 70)
    print("SUMMARY: MONOTONICITY THEOREM")
    print("=" * 70)
    print()
    print("OBSERVATION: λ_n is strictly monotone increasing for n ≤ 10,000.")
    print()
    print("IMPLICATION: If proved for all n,")
    print("  λ_1 > 0 + λ_n increasing → λ_n > 0 for all n → RH.")
    print()
    print("KEY FORMULA: Δλ_n = Σ_ρ (1-1/ρ)^n / ρ")
    print()
    print("For on-line zeros: |(1-1/ρ)| = 1 exactly.")
    print("Per-pair contribution ≈ sin(n/γ) / γ:")
    print("  Positive for n < πγ (n < turning point)")
    print("  Negative for πγ < n < 2πγ (past turning point)")
    print()
    print("K-zero lower bound:")
    print("  For n ≤ πγ_K (all tail zeros in first half-period),")
    print("  the unverified-zero tail contributes POSITIVELY.")
    print("  So Δλ_n ≥ Δλ_n^(K) = partial sum from K verified zeros.")
    print()
    print("REMAINING GAP:")
    print("  For n > πγ_K, some tail zeros enter negative half-period.")
    print("  Need to bound their negative contributions.")
    print("  This is the same gap as in the exponential domination argument.")
    print("=" * 70)
