"""
Weil Explicit Formula — Sign Analysis for λ_n

The Weil explicit formula gives:
  λ_n = Σ_ρ [1 - (1 - 1/ρ)^n]

where the sum is over ALL nontrivial zeros ρ = σ + iγ.

Each quadruple {ρ, ρ̄, 1-ρ, 1-ρ̄} contributes to λ_n.
By the compact pair formula, the quadruple contribution is:
  Q_n(σ,γ) = F_n(σ,γ) = 4(1 - cosh(n·log r₁(σ)) · cos(nθ(σ,γ)))

Key structural fact:
  - Q_n(1/2, γ) = 4(1 - cos(nθ(1/2,γ))) ≥ 0  [always safe on line]
  - Q_n(σ, γ) < 0 requires σ ≠ 1/2 AND resonance nθ ≈ 2πk

The question: can an off-line zero at σ₀ + iγ₀ (with σ₀ ≠ 1/2) cause λ_n < 0
for some n, given all other zeros are on the line?

SINGLE OFF-LINE ZERO SCENARIO:
  Replace one zero ρ* = σ₀ + iγ₀ (and its conjugates/reflections) with
  an on-line zero at 1/2 + iγ₀. The change in λ_n is:
    ΔQ_n = F_n(σ₀, γ₀) - F_n(1/2, γ₀)

  If ΔQ_n < 0 for some n, then moving a zero off-line DECREASES λ_n.
  We want to know: can the decrease be large enough to make λ_n < 0?

ANALYSIS:
  F_n(1/2, γ) = 4(1 - cos(nφ)) where φ = π - 2·arctan(2γ) ≈ 1/γ for large γ
  F_n(σ, γ) = 4(1 - cosh(n·log r₁(σ)) · cos(nθ(σ,γ)))

  The "damage" from moving ρ* off-line at resonance is approximately:
    ΔQ_n ≈ 4(1 - C · (-1)) - 4(1 - 1·(+1)) = 4(1+C) - 4(0) = 4(1+C)
  Wait — at resonance with a SAFE on-line zero:
    On-line: cos(nφ) ≈ 1 → F_n(1/2,γ) ≈ 0 (near resonance for on-line zero)
    Off-line at resonance: cosh >> 1, cos(nθ) > 0 → F_n(σ,γ) < 0

  So moving a zero off-line at a resonance causes:
    ΔQ_n ≈ F_n(σ,γ) - F_n(1/2,γ)  (negative at resonance)

This script:
1. Computes ΔQ_n for a hypothetical off-line zero at various (σ,γ,n)
2. Finds the maximum negative contribution from a single off-line zero
3. Compares to the sum of positive contributions from all OTHER on-line zeros
4. Asks: can we prove the sum is always positive?
"""

import numpy as np

GAMMA_1 = 14.134725141734693790


def theta(sigma, gamma):
    """θ(σ,γ) = arctan((1-σ)/γ) + arctan(σ/γ)."""
    return np.arctan((1 - sigma) / gamma) + np.arctan(sigma / gamma)


def log_r1(sigma, gamma):
    """log|α₁(σ,γ)| = (1/2)·log(((1-σ)²+γ²)/(σ²+γ²))."""
    return 0.5 * np.log(((1 - sigma)**2 + gamma**2) / (sigma**2 + gamma**2))


def F_pair(n, sigma, gamma):
    """F_n(σ,γ) = 4(1 - cosh(n·log r₁) · cos(nθ)) — compact formula."""
    lr = log_r1(sigma, gamma)
    th = theta(sigma, gamma)
    return 4.0 * (1.0 - np.cosh(n * lr) * np.cos(n * th))


def F_online(n, gamma):
    """F_n(1/2, γ) = 4(1 - cos(n·φ)), φ = θ(1/2,γ) = π - 2·arctan(2γ)."""
    phi = theta(0.5, gamma)
    return 4.0 * (1.0 - np.cos(n * phi))


def damage_from_off_line(n, sigma, gamma):
    """
    ΔQ_n = F_n(σ,γ) - F_n(1/2,γ).
    Negative means off-line zero is WORSE than on-line at this n.
    """
    return F_pair(n, sigma, gamma) - F_online(n, gamma)


def find_worst_resonance(sigma, gamma, n_max=10000):
    """
    Find the n ≤ n_max where a zero at (σ,γ) gives the most negative Q_n.
    Returns (worst_n, worst_Qn, worst_ΔQn).
    """
    worst_n = None
    worst_Q = 0.0
    worst_dQ = 0.0
    for n in range(1, n_max + 1):
        Q = F_pair(n, sigma, gamma)
        dQ = damage_from_off_line(n, sigma, gamma)
        if Q < worst_Q:
            worst_Q = Q
            worst_n = n
            worst_dQ = dQ
    return worst_n, worst_Q, worst_dQ


def on_line_contribution_sum(n, gammas):
    """Sum of F_n(1/2,γ) over all on-line zeros (positive by definition)."""
    total = 0.0
    for g in gammas:
        total += F_online(n, g)
    return total


def analyze_perturbation(sigma_off, gamma_off, gammas_online, n_range=(1, 5000)):
    """
    Suppose all zeros are on the line EXCEPT one at (sigma_off, gamma_off).
    What is the minimum of:
      Σ_{j≠*} F_n(1/2,γ_j) + F_n(sigma_off, gamma_off)
    over n?
    """
    n_vals = np.arange(n_range[0], n_range[1] + 1)

    # Online contributions from all-but-one zeros
    online_sums = np.array([on_line_contribution_sum(n, gammas_online) for n in n_vals])

    # Off-line contribution from the perturbed zero
    offline_contr = np.array([F_pair(n, sigma_off, gamma_off) for n in n_vals])

    total = online_sums + offline_contr
    min_idx = np.argmin(total)
    return n_vals[min_idx], total[min_idx], online_sums[min_idx], offline_contr[min_idx]


def resonance_map(sigma, gamma, n_max=1000):
    """Show F_n(σ,γ) and F_n(1/2,γ) side by side for resonance analysis."""
    print(f"\nResonance map for σ={sigma:.3f}, γ={gamma:.4f}:")
    print(f"{'n':>6} {'Q_off':>12} {'Q_on':>12} {'ΔQ':>12} {'note':>10}")
    phi = theta(0.5, gamma)
    th0 = theta(sigma, gamma)
    for n in range(1, n_max + 1):
        Q_off = F_pair(n, sigma, gamma)
        Q_on = F_online(n, gamma)
        dQ = Q_off - Q_on
        # Flag resonances: nθ near 2πk
        residue = (n * th0) % (2 * np.pi)
        near_res = residue < 0.3 or residue > 2 * np.pi - 0.3
        if Q_off < -0.1 or near_res:
            k = round(n * th0 / (2 * np.pi))
            print(f"{n:>6} {Q_off:>12.4f} {Q_on:>12.4f} {dQ:>12.4f}  ← nθ≈{k}×2π")


def size_of_damage(sigma_range=(0.01, 0.49), gamma=14.1347, n=89):
    """
    For a given (n, γ), how bad can F_n(σ,γ) get as σ varies?
    The cosh blows up as σ→0 (r₁→r₁_max).
    """
    sigmas = np.linspace(sigma_range[0], sigma_range[1], 500)
    Fvals = np.array([F_pair(n, s, gamma) for s in sigmas])
    print(f"\nF_{n}(σ, γ={gamma:.4f}) as σ varies:")
    print(f"  Min: {Fvals.min():.4f} at σ={sigmas[Fvals.argmin()]:.4f}")
    print(f"  Max: {Fvals.max():.4f}")
    print(f"  F at σ=1/2: {F_online(n, gamma):.4f}")
    print(f"  cosh(n·log r₁(0,γ)): {np.cosh(n * 0.5 * np.log(1 + 1/gamma**2)):.4f}")
    print(f"  cos(nθ(0,γ)): {np.cos(n * theta(0, gamma)):.6f}")


def structural_inequality_test(n, gammas, sigma_off, gamma_off):
    """
    TEST: λ_n > 0 even with one off-line zero.
    Lower bound:
      λ_n ≥ Σ_j F_n(1/2,γ_j) + F_n(σ_off,γ_off)
            [using F_n(1/2) as lower bound for on-line zeros]
    """
    online_sum = sum(F_online(n, g) for g in gammas)
    offline = F_pair(n, sigma_off, gamma_off)
    lb = online_sum + offline
    print(f"\nStructural inequality test: n={n}, off-line at (σ={sigma_off:.3f}, γ={gamma_off:.4f})")
    print(f"  Online sum (all {len(gammas)} zeros at 1/2): {online_sum:.4f}")
    print(f"  Off-line contribution:                       {offline:.4f}")
    print(f"  Lower bound on λ_n:                         {lb:.4f}")
    print(f"  Result: {'λ_n > 0 ✓' if lb > 0 else 'NEGATIVE — possible violation ✗'}")
    return lb


if __name__ == "__main__":
    print("=" * 70)
    print("Weil Explicit Formula — Sign Analysis for λ_n")
    print("=" * 70)

    # The most dangerous scenario: n=89 is near the first resonance of γ₁
    # (since 89 × arctan(1/14.135) ≈ 89 × 0.07063 ≈ 6.286 ≈ 2π)
    n_res = 89
    sigma_off = 0.05   # hypothetical off-line zero close to real axis
    gamma_off = GAMMA_1

    print(f"\n--- WORST CASE: n={n_res}, off-line at σ={sigma_off}, γ=γ₁ ---")
    print(f"nθ(σ=0.05, γ₁) = {n_res * theta(sigma_off, gamma_off):.4f}")
    print(f"2π = {2*np.pi:.4f}")
    print(f"cosh(n·log r₁) = {np.cosh(n_res * log_r1(sigma_off, gamma_off)):.4f}")
    print(f"F_{n_res}(0.05, γ₁) = {F_pair(n_res, sigma_off, gamma_off):.4f}")
    print(f"F_{n_res}(1/2, γ₁)  = {F_online(n_res, gamma_off):.4f}")

    # Size of damage for various σ at resonance n
    size_of_damage(gamma=GAMMA_1, n=n_res)

    # How many on-line zeros needed to compensate worst off-line damage?
    worst_offline = min(F_pair(n_res, s, GAMMA_1)
                        for s in np.linspace(0.001, 0.499, 1000))
    print(f"\nWorst F_{n_res}(σ,γ₁) over σ∈[0,1/2]: {worst_offline:.4f}")
    print(f"To compensate, need: ≥ {abs(worst_offline) / F_online(n_res, GAMMA_1):.1f} × F(1/2,γ₁)")
    print(f"  (there are infinitely many other zeros contributing positively)")

    # Key asymptotic analysis
    print("\n--- ASYMPTOTIC GROWTH ANALYSIS ---")
    print("As n → ∞, how does the sum of online contributions grow?")
    print("Online: Σ_γ 4(1 - cos(nφ_γ)) → 4N_zeros × 1 = 4N_zeros (on average)")
    print("The sum grows as ~4N_zeros, while worst offline ≈ -4·cosh(n·log r₁)")
    print("cosh(n·log r₁) grows EXPONENTIALLY with n for σ ≠ 1/2")
    print("But the DENSITY of zeros N(T) ~ T·log(T)/(2π) grows too")
    print()

    # Demonstrate the exponential vs linear growth tension
    print("Exponential (cosh) vs linear (zero count) growth:")
    print(f"{'n':>8} {'cosh(n·logr₁)':>18} {'~4n zeros up to n/γ₁':>22}")
    for n in [10, 50, 100, 500, 1000, 3000, 5000, 10000]:
        # r₁ at σ=0.1, γ=γ₁ (a reasonably off-line zero)
        lr = log_r1(0.1, GAMMA_1)
        c = np.cosh(n * lr)
        # rough count of zeros γ ≤ n (since φ = θ(1/2,γ) ≈ 1/γ, resonance at n/γ ≈ 2πk)
        approx_zeros = n * np.log(n) / (2 * np.pi) if n > 1 else 1
        print(f"{n:>8} {c:>18.2e} {approx_zeros:>22.0f}")

    print("\nConclusion: cosh grows as e^{n·|log r₁|}, while number of zeros")
    print("grows as n·log(n). For any fixed σ≠1/2, cosh WINS for large n.")
    print("This means a single off-line zero can make λ_n < 0 for large enough n.")
    print("→ This confirms RH structurally: off-line zeros ARE incompatible with λ_n > 0.")
    print()
    print("The challenge: proving this formally = proving no off-line zeros exist.")
    print("Circular unless we find a fundamentally different angle.")

    print("\n--- WHAT IF WE BOUND THE GROWTH RATE? ---")
    print("Suppose σ* ∉ [1/2-ε, 1/2+ε] for some ε > 0.")
    print("Then |log r₁(σ*)| ≥ δ(ε) > 0 (uniform lower bound).")
    print("Then F_n(σ*,γ*) ≤ 4(1 - e^{nδ}/2) < -e^{nδ} for large n.")
    print("But Σ_{on-line} F_n(1/2,γ) ≤ Σ 4·2 = 8·N_zeros (each term ≤ 8)")
    print("→ For n large enough: ONE off-line zero overwhelms ALL on-line zeros.")
    print("→ λ_n < 0 for large n if ANY zero is off-line by ε.")
    print()
    print("This IS the proof structure — but we need:")
    print("  1. A lower bound on δ(ε) = min|log r₁| for σ off-line by ε")
    print("  2. An upper bound on the sum Σ F_n(1/2,γ) that doesn't grow exponentially")
    print("  3. The n threshold above which the off-line term dominates")
    print()
    print("Alternatively: use the ACTUAL λ_n formula (not just lower bounds)")
    print("and exploit that λ_n > 0 for ALL n ≤ 4.83×10¹² (Platt-Trudgian).")

    print("\n--- EXPONENTIAL DOMINATION THEOREM (sketch) ---")
    # If σ* is off-line at distance ε from 1/2, quantify when cosh dominates.
    eps = 0.1  # off-line by 0.1
    sigma_star = 0.5 - eps  # = 0.4
    gamma_star = GAMMA_1
    log_r1_val = log_r1(sigma_star, gamma_star)
    print(f"Off-line zero at σ=0.4 (ε=0.1 from 1/2), γ=γ₁:")
    print(f"  log r₁(0.4, γ₁) = {log_r1_val:.6f}")
    print(f"  So cosh(n·log r₁) ≥ (1/2)·e^{{n·{log_r1_val:.4f}}}")
    print(f"  Worst F_n term ≈ -4·(1/2)·e^{{nδ}} where δ = {log_r1_val:.4f}")
    print()
    print(f"  Online sum upper bound: Σ_γ 4·2 = 8N(T) where N(T) ~ T·log(T)/(2π)")
    print(f"  For the off-line term to dominate: e^{{nδ}} >> N(T)")
    print(f"  But T is the height of the off-line zero = γ* ≈ {GAMMA_1:.2f}")
    print(f"  N({GAMMA_1:.2f}) ≈ {int(GAMMA_1 * np.log(GAMMA_1) / (2*np.pi))} zeros")
    print()
    threshold_n = np.log(1000) / log_r1_val  # when cosh > 1000 >> N(γ*)
    print(f"  cosh(nδ) > 1000 when n > {threshold_n:.1f}")
    print(f"  At n={int(threshold_n)+1}: cosh = {np.cosh((int(threshold_n)+1)*log_r1_val):.2e}")
    print(f"  At n={int(threshold_n)+1}: online sum ≤ {8 * int(GAMMA_1 * np.log(GAMMA_1) / (2*np.pi)):.0f}")
    print()
    print("  PROBLEM: this only counts zeros up to γ*. But λ_n sums over ALL zeros.")
    print("  The full online sum Σ_ALL F_n(1/2,γ) diverges — it needs regularization.")
    print("  The Weil formula's sum is conditionally convergent, not absolutely convergent.")
    print("  This is the real analytic subtlety.")
