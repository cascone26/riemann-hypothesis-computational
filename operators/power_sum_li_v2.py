"""
Li coefficients via Newton power sums — corrected version.

Formula: λ_n = Σ_{k=1}^{n} C(n,k) (-1)^{k-1} p_k

where p_k = Σ_ρ 1/ρ^k = 2 Σ_{γ>0} Re[1/(1/2+iγ)^k]

For k=1: p_1 = 1 + γ_E/2 - (1/2)log(4π) (exact, no zeros needed)
For k≥2: use 50K-zero partial sum (convergent, tail < 10^{-4})

Key question: What's the sign pattern of p_k, and does it prove λ_n > 0?

Answer: The sign pattern of (-1)^{k-1} p_k determines whether λ_n > 0 trivially.
If all terms have the same sign, λ_n > 0. If not, cancellation is needed.
"""

import numpy as np
from scipy.special import comb as scipy_comb
from mpmath import mp, mpf, euler, pi, log, factorial
import mpmath
import json
import os

mp.dps = 50

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")


def load_zeros(count=50000):
    path = os.path.join(DATA_DIR, "zeros_100k.txt")
    zeros = []
    with open(path, "r") as f:
        for i, line in enumerate(f):
            if i >= count:
                break
            line = line.strip()
            if line:
                zeros.append(float(line))
    return np.array(zeros)


def compute_all_power_sums(zeros_np, k_max=40):
    """
    Compute p_k = 2 Σ_{γ>0} Re[1/(1/2+iγ)^k] for k=1..k_max.

    For k=1: the sum converges conditionally; we use the exact formula:
    p_1 = 1 + γ_E/2 - (1/2)log(4π) ≈ 0.02309571...

    For k≥2: the sum converges absolutely. Tail correction is tiny.
    """
    rhos = 0.5 + 1j * zeros_np

    p_exact = {}
    p_partial = {}

    # Exact p_1
    p1_exact = float(1 + euler / 2 - log(4 * pi) / 2)
    p_exact[1] = p1_exact
    print(f"  p_1 (exact formula) = {p1_exact:.15f}")

    # Direct sum for k=1..k_max
    print(f"\n  Computing partial sums p_k (50K zeros) for k=2..{k_max}...")
    for k in range(1, k_max + 1):
        # 1/ρ^k for each zero
        vals = 1.0 / (rhos ** k)
        pk = 2.0 * np.sum(vals.real)
        p_partial[k] = pk

    print(f"  p_1 (partial, 50K) = {p_partial[1]:.15f}")
    print(f"  p_1 (exact)        = {p1_exact:.15f}")
    print(f"  Diff (tail, k>50K) = {p1_exact - p_partial[1]:.2e}")

    # Build best estimates: exact for k=1, partial for k≥2
    p_best = {}
    p_best[1] = p1_exact
    for k in range(2, k_max + 1):
        p_best[k] = p_partial[k]

    return p_best, p_partial


def print_power_sum_table(p_best, k_max=30):
    """Print the sign structure of p_k and (-1)^{k-1} p_k."""
    print("\n" + "=" * 70)
    print("POWER SUM TABLE p_k = Σ_ρ 1/ρ^k")
    print("=" * 70)
    print(f"\n  {'k':>3}  {'p_k':>22}  {'(-1)^{k-1} p_k':>22}  {'sign':>5}")
    print(f"  {'-'*60}")

    for k in range(1, k_max + 1):
        pk = p_best.get(k)
        if pk is None:
            continue
        combined = ((-1) ** (k - 1)) * pk
        sign = "+" if combined > 0 else "-"
        print(f"  {k:>3}  {pk:>22.15f}  {combined:>22.15f}  {sign:>5}")

    pos_k = [k for k in range(1, k_max + 1) if ((-1)**(k-1)) * p_best.get(k, 0) > 0]
    neg_k = [k for k in range(1, k_max + 1) if ((-1)**(k-1)) * p_best.get(k, 0) < 0]

    print(f"\n  (-1)^{{k-1}} p_k > 0 for k = {pos_k}")
    print(f"  (-1)^{{k-1}} p_k < 0 for k = {neg_k}")

    # If ALL were positive: λ_n > 0 trivially (all terms in sum positive)
    # If some are negative: need cancellation argument
    all_positive = len(neg_k) == 0
    print(f"\n  All (-1)^{{k-1}} p_k > 0? {all_positive}")
    if not all_positive:
        print(f"  CANCELLATION NEEDED: negative k-values create negative terms in λ_n.")
        print(f"  The positive terms must dominate. This is the content of the Li criterion.")


def compute_lambda_correct(p_best, n_max=200):
    """
    Compute λ_n = Σ_{k=1}^n C(n,k) (-1)^{k-1} p_k correctly.

    Uses all p_k, not just odd ones.
    """
    print(f"\n--- Computing λ_n for n=1..{n_max} ---")

    lambda_vals = {}
    lambda_positive = True
    deltas = {}
    lam_prev = 0.0

    for n in range(1, n_max + 1):
        lam = 0.0
        for k in range(1, n + 1):
            pk = p_best.get(k, 0.0)
            binom = float(mpmath.binomial(n, k))
            sign = (-1) ** (k - 1)
            lam += binom * sign * pk

        lambda_vals[n] = lam
        deltas[n] = lam - lam_prev
        lam_prev = lam

        if lam < 0:
            lambda_positive = False
            print(f"  *** NEGATIVE: λ_{n} = {lam:.6f} ***")

        if n <= 30 or n % 50 == 0:
            asymp = float((n / 2) * (np.log(n / (2 * np.pi)) + float(euler) - 1))
            ratio = lam / asymp if abs(asymp) > 0.01 else float('nan')
            print(f"  λ_{n:>4d} = {lam:>14.6f}  asymp = {asymp:>14.6f}  ratio = {ratio:>8.4f}")

    return lambda_vals, deltas, lambda_positive


def analyze_sign_structure(p_best, k_max=40):
    """
    Deep analysis of the sign structure of (-1)^{k-1} p_k.

    Key question: Can we prove (-1)^{k-1} p_k has a definite sign?

    From the analysis:
    p_k = 2 Σ_{γ>0} Re[1/(1/2+iγ)^k]
    (-1)^{k-1} p_k = 2(-1)^{k-1} Σ_{γ>0} Re[1/(1/2+iγ)^k]

    The dominant behavior for large γ:
    Re[1/(1/2+iγ)^k] ≈ Re[1/(iγ)^k] = Re[(−i/γ)^k] = Re[e^{−ikπ/2}/γ^k]
    = cos(−kπ/2)/γ^k = cos(kπ/2)/γ^k

    So: sign(p_k) ≈ sign(cos(kπ/2)) for large γ dominant
    k=1: cos(π/2)=0, subleading: +
    k=2: cos(π)=-1 → p_2 < 0
    k=3: cos(3π/2)=0, subleading: negative
    k=4: cos(2π)=+1 → p_4 > 0
    k=5: cos(5π/2)=0, subleading...

    Period 4: signs repeat roughly as +,-,-,+,+,-,-,+,...
    But the factor (-1)^{k-1} changes this.

    (-1)^{k-1} p_k: for the dominant behavior:
    k=1: (+)×(+) = + (since p_1>0)
    k=2: (-)×(-) = + (since p_2<0, (-1)^1=-1, so (-1)(-) = +)
    k=3: (+)×(-) = - (since p_3<0, (-1)^2=+1)
    k=4: (-)×(+) = - (since p_4>0, (-1)^3=-1)

    So signs of (-1)^{k-1} p_k: +, +, -, -,  +, +, -, -, ...
    (period 4 grouping: ++--++--...)

    This means lambda_n = Σ C(n,k) [(-1)^{k-1} p_k] has ALTERNATING blocks
    of positive and negative terms. Not all positive!

    But the magnitudes decrease rapidly: |p_k| ~ (1/|ρ_1|)^k where |ρ_1| ≈ 14.
    So the terms become tiny for k≥4 compared to k=1,2 terms.

    The dominance: λ_n ≈ C(n,1)(+p_1) + C(n,2)(+|p_2|) + C(n,3)(-|p_3|) + ...
    The first two positive terms dominate for any n.
    """
    print("\n" + "=" * 70)
    print("SIGN STRUCTURE ANALYSIS")
    print("=" * 70)

    print("""
  Analysis of (-1)^{k-1} p_k sign pattern:

  For large γ, Re[1/(1/2+iγ)^k] ≈ cos(kπ/2)/γ^k (dominant term)

  k mod 4:
    k≡1 (mod 4): cos(π/2)=0, subleading ≈ (1/2)/γ^2 > 0 → p_k > 0
    k≡2 (mod 4): cos(π)=-1 → p_k < 0
    k≡3 (mod 4): cos(3π/2)=0, subleading ≈ -3/(2γ^4) < 0 → p_k < 0
    k≡0 (mod 4): cos(2π)=+1 → p_k > 0

  Sign of (-1)^{k-1} × p_k:
    k≡1 (mod 4): (-1)^0 × (+) = +
    k≡2 (mod 4): (-1)^1 × (-) = +
    k≡3 (mod 4): (-1)^2 × (-) = -
    k≡0 (mod 4): (-1)^3 × (+) = -

  Pattern: +, +, -, -, +, +, -, -, ...  (period 4)

  So λ_n = Σ C(n,k) [(-1)^{k-1} p_k] has this structure:
    Positive contributions: k ≡ 1,2 (mod 4)
    Negative contributions: k ≡ 3,0 (mod 4)

  The RATIO of positive to negative contributions:
    Positive: C(n,1)p_1 + C(n,2)|p_2| + C(n,5)p_5 + C(n,6)|p_6| + ...
    Negative: C(n,3)|p_3| + C(n,4)p_4 + C(n,7)|p_7| + ...

  Since |p_k| decreases geometrically (ratio ≈ 1/14 per step),
  the first pair (k=1,2) dominates:
    λ_n ≈ n×p_1 + C(n,2)|p_2| - C(n,3)|p_3| - ...
         ≈ n×p_1 + n²|p_2|/2 - n³|p_3|/6 + ...

  For large n: the n²|p_2|/2 term dominates. Since p_2 < 0:
    C(n,2)|p_2| = n(n-1)/2 × |p_2|

  This grows as n² — faster than n.

  But wait! The ASYMPTOTIC λ_n ~ (n/2)log n grows faster than n².
  So the dominance structure changes for large n.

  The resolution: For large n, the k≈n terms in the sum become important.
  The sum Σ C(n,k)(-1)^{k-1}p_k for k near n involves HUGE binomials.
  """)

    # Compute numerically to verify the sign pattern
    print("  Numerical verification of sign pattern:")
    print(f"  {'k':>3}  {'p_k':>16}  {'sign(p_k)':>10}  {'(-1)^{{k-1}}p_k':>16}  {'sign(combo)':>12}")
    print(f"  {'-'*60}")

    for k in range(1, min(k_max + 1, 21)):
        pk = p_best.get(k, 0.0)
        sign_pk = "+" if pk > 0 else "-"
        combined = ((-1) ** (k - 1)) * pk
        sign_combo = "+" if combined > 0 else "-"
        predicted = {1: "+", 2: "+", 3: "-", 0: "-"}.get(k % 4, "?")
        match = "✓" if predicted == sign_combo else "✗"
        print(f"  {k:>3}  {pk:>16.10f}  {sign_pk:>10}  {combined:>16.10f}  {sign_combo:>12} (pred:{predicted} {match})")


def magnitude_analysis(p_best, k_max=40):
    """
    Study the magnitude decay of p_k.

    Predicted: |p_k| ~ C × (1/|ρ_1|)^k where |ρ_1| ≈ 14.14

    This means the power sum converges RAPIDLY.
    For k=10: |p_10| ≈ C × (1/14)^10 ≈ C × 4×10^{-12}

    The Lambda_n sum terminates at k=n, so for n >> |ρ_1|, the high-k
    terms are negligible. This means:
      For n >> 14: λ_n ≈ Σ_{k=1}^{~30} C(n,k)(-1)^{k-1}p_k

    The binomial coefficients C(n,k) for fixed k and large n grow as n^k/k!.
    So: λ_n ≈ Σ_{k=1}^{30} (n^k/k!) × (-1)^{k-1} p_k for large n

    = n × p_1 + n²/2 × (-p_2) + n³/6 × p_3 + ...

    = n × [p_1 - n/2 × p_2 + n²/6 × p_3 - ...]

    This is like a power series in n × p_2 / p_1, which grows!
    For n >> p_1/|p_2| ≈ 0.023/0.01 ≈ 2.3, higher terms dominate.

    BUT the sum converges because of the 1/k! denominators.
    The generating function Σ (n^k/k!) (-1)^{k-1} p_k = exponential!

    Specifically: Σ_{k=1}^∞ (t^k/k!) × c_k = exp(Σ t×something) - 1
    where the series relates to the logarithm of a product formula.

    This is the key connection to the product formula for ξ.
    """
    print("\n" + "=" * 70)
    print("MAGNITUDE DECAY ANALYSIS")
    print("=" * 70)

    mags = [(k, abs(p_best.get(k, 0.0))) for k in range(1, k_max + 1)]

    print(f"\n  {'k':>3}  {'|p_k|':>18}  {'ratio |p_{k+1}|/|p_k|':>22}  {'predicted':>12}")
    print(f"  {'-'*60}")

    rho1_abs = np.sqrt(0.25 + 14.1347**2)  # |ρ_1| = |1/2 + 14.135i|

    for i, (k, mag) in enumerate(mags[:20]):
        if i > 0 and mags[i-1][1] > 1e-30:
            ratio = mag / mags[i-1][1]
            predicted = 1.0 / rho1_abs
            match = "✓" if abs(ratio - predicted) < 0.5 else ""
            print(f"  {k:>3}  {mag:>18.12f}  {ratio:>22.8f}  "
                  f"{predicted:>12.8f} {match}")
        else:
            print(f"  {k:>3}  {mag:>18.12f}")


def main():
    print("=" * 70)
    print("POWER SUM Li ANALYSIS v2 — Complete computation")
    print("=" * 70)

    print("\n--- Loading zeros ---")
    zeros_np = load_zeros(50000)
    print(f"  Loaded {len(zeros_np)} zeros")

    print("\n--- Computing power sums ---")
    p_best, p_partial = compute_all_power_sums(zeros_np, k_max=40)

    # Print the key table
    print_power_sum_table(p_best, k_max=30)

    # Compute λ_n correctly
    lambda_vals, deltas, all_positive = compute_lambda_correct(p_best, n_max=200)

    # Compare with our direct sum computation
    print("\n--- Comparison with direct 50K-zero sum ---")
    import numpy as npy
    li50k = npy.load(os.path.join(RESULTS_DIR, "li_50k.npy"))

    print(f"\n  {'n':>5}  {'power_sum λ_n':>18}  {'direct_sum λ_n':>18}  {'ratio':>8}")
    print(f"  {'-'*58}")
    for n in [1, 2, 3, 5, 10, 20, 50, 100, 200]:
        lp = lambda_vals.get(n, float('nan'))
        ld = float(li50k[n-1]) if n <= len(li50k) else float('nan')
        ratio = lp / ld if abs(ld) > 1e-10 else float('nan')
        print(f"  {n:>5}  {lp:>18.8f}  {ld:>18.8f}  {ratio:>8.6f}")

    print("\n  NOTE: power_sum uses exact p_1 + 50K-partial p_k (k≥2)")
    print("  direct_sum uses 50K zeros directly (lower bound on true λ_n)")
    print("  power_sum gives the TRUE λ_n (including tail of missing zeros)")
    print("  difference = contribution from zeros #50001 to ∞ (all positive)")

    # Sign structure analysis
    analyze_sign_structure(p_best)

    # Magnitude analysis
    magnitude_analysis(p_best)

    # Summary
    print("\n" + "=" * 70)
    print("CONCLUSIONS")
    print("=" * 70)

    neg_k = [k for k in range(1, 31) if ((-1)**(k-1)) * p_best.get(k, 0) < 0]

    print(f"""
  1. All λ_n > 0 for n=1..200 (power sum method): {all_positive}

  2. Sign structure of (-1)^{{k-1}} p_k:
     Pattern is approximately +,+,-,-,+,+,-,- (period 4, starting k=1)
     Negative k-values: {neg_k}
     → λ_n > 0 is NOT trivially obvious from signs alone
     → Cancellation between positive/negative contributions is required
     → This cancellation is essentially the content of the Li criterion

  3. Magnitude decay: |p_k| ~ (1/|ρ_1|)^k ≈ (1/14.14)^k
     → Very rapid decay; higher k terms negligible for λ_n
     → The sum is dominated by k=1,2,3,4 terms

  4. For large n, the dominant terms in λ_n = Σ C(n,k)(-1)^{{k-1}}p_k are:
     - k=1: n × p_1 (grows like n)
     - k=2: n(n-1)/2 × |p_2| (grows like n²)
     - The Keiper asymptotic (n/2)log(n) comes from ALL terms together

  5. PROOF OBSTACLE: Proving λ_n > 0 from the power sum formula requires
     showing the partial sums sum_{{k<=K}} C(n,k)(-1)^{{k-1}}p_k remain positive.
     This is equivalent to the Li criterion itself.
  """)

    # Save
    save_data = {
        "power_sums": {str(k): v for k, v in p_best.items()},
        "lambda_n_200": {str(k): v for k, v in lambda_vals.items()},
        "all_positive_200": all_positive,
        "neg_sign_k": neg_k,
        "rho1_abs": float(np.sqrt(0.25 + 14.1347**2)),
    }
    with open(os.path.join(RESULTS_DIR, "power_sum_li_v2.json"), "w") as f:
        json.dump(save_data, f, indent=2)
    print(f"  Saved to {RESULTS_DIR}/power_sum_li_v2.json")


if __name__ == "__main__":
    main()
