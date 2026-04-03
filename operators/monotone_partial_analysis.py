"""
Monotonicity and partial sum stability analysis of the Li criterion.

1. MONOTONICITY: Is lambda_{n+1} - lambda_n always positive?
   If so, lambda_1 = p_1 > 0 (unconditional) is the global minimum.

2. PARTIAL SUM STABILITY: Does S_n(K) = sum_{k=1}^K C(n,k)(-1)^{k-1} p_k
   stay positive for all K <= n?
   If yes: even with only the first K power sums, lambda_n > 0.

3. PHASE FORMULA: lambda_n = 2 * sum_gamma [1 - cos(n * phi_gamma)]
   where phi_gamma = pi - 2*arctan(2*gamma) > 0.
   Each term >= 0, so lambda_n >= 0 trivially under RH.
   The question: does this positivity follow from the phase structure alone?
"""

import numpy as np
import json
import os
import sys
from scipy.special import comb as scipy_comb

RESULTS_DIR = "/Users/scones/projects/riemann/results"
ZEROS_FILE = "/Users/scones/projects/riemann/data/zeros_odlyzko_100k"

def load_power_sums():
    """Load p_k values from saved results."""
    path = os.path.join(RESULTS_DIR, "power_sum_li_v2.json")
    with open(path) as f:
        data = json.load(f)
    p = {}
    for k_str, v in data["power_sums"].items():
        p[int(k_str)] = v
    return p

def load_li_array():
    """Load the full 50K Li array."""
    path = os.path.join(RESULTS_DIR, "li_50k.npy")
    return np.load(path)

def load_zeros(n=50000):
    """Load zeros from text file."""
    path = ZEROS_FILE
    # Try loading as text
    data = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    data.append(float(line.split()[0]))
                except:
                    pass
            if len(data) >= n:
                break
    return np.array(data)

def phase_from_gamma(gamma):
    """phi_gamma = pi - 2*arctan(2*gamma). Always in (0, pi) for gamma > 0."""
    return np.pi - 2.0 * np.arctan(2.0 * gamma)

def main():
    print("=" * 70)
    print("MONOTONICITY AND PARTIAL SUM STABILITY ANALYSIS")
    print("=" * 70)

    # =========================================================================
    # 1. MONOTONICITY: Check delta_lambda_n = lambda_{n+1} - lambda_n
    # =========================================================================
    print("\n--- Monotonicity of lambda_n (from 50K stored array) ---")
    li = load_li_array()  # li[k] = lambda_{k+1} (0-indexed)
    n_max = len(li)

    deltas = np.diff(li)  # deltas[k] = li[k+1] - li[k] = lambda_{k+2} - lambda_{k+1}

    neg_delta_count = np.sum(deltas < 0)
    min_delta = np.min(deltas)
    min_delta_n = np.argmin(deltas) + 1  # n where delta is smallest (lambda_{n+1} - lambda_n)

    print(f"  Total deltas computed: {len(deltas)}")
    print(f"  Number of negative deltas: {neg_delta_count}")
    print(f"  Min delta: {min_delta:.6f} at n={min_delta_n}")
    print(f"  Max delta: {np.max(deltas):.4f}")
    print(f"  Monotone increasing? {neg_delta_count == 0}")

    # Show small deltas
    print(f"\n  First 20 deltas (lambda_{{n+1}} - lambda_n):")
    print(f"  {'n':>6}  {'delta':>12}  {'lambda_n':>12}")
    for n in range(1, 21):
        if n-1 < len(deltas):
            print(f"  {n:>6}  {deltas[n-1]:>12.6f}  {li[n-1]:>12.6f}")

    # =========================================================================
    # 2. PARTIAL SUM STABILITY using p_k values
    # =========================================================================
    print("\n--- Partial Sum Stability S_n(K) ---")
    print("  (S_n(K) = sum_{k=1}^K C(n,k)*(-1)^{k-1}*p_k)")

    p = load_power_sums()
    k_max = max(p.keys())
    print(f"  Available p_k: k=1..{k_max}")

    # Check stability for n=1..30 (limited by available p_k)
    n_test = min(30, k_max)

    min_partial = float('inf')
    min_partial_n = -1
    min_partial_K = -1
    dip_count = 0

    print(f"\n  Checking n=1..{n_test}:")
    print(f"  {'n':>4}  {'K_dip':>6}  {'min_partial':>14}  {'lambda_n':>12}  {'monotone':>8}")

    for n in range(1, n_test + 1):
        partial_sums = []
        running = 0.0
        for k in range(1, n + 1):
            pk = p.get(k, 0.0)
            sign = (-1) ** (k - 1)
            binom = scipy_comb(n, k, exact=True)
            running += binom * sign * pk
            partial_sums.append(running)

        min_ps = min(partial_sums)
        min_ps_K = partial_sums.index(min_ps) + 1
        final = partial_sums[-1]

        # Is it always positive?
        monotone = all(ps > 0 for ps in partial_sums)

        if min_ps < min_partial:
            min_partial = min_ps
            min_partial_n = n
            min_partial_K = min_ps_K

        if min_ps < 0:
            dip_count += 1

        status = "YES" if monotone else "NO"
        print(f"  {n:>4}  {min_ps_K:>6}  {min_ps:>14.8f}  {final:>12.8f}  {status:>8}")

    print(f"\n  Summary:")
    print(f"  Cases where S_n(K) dipped negative: {dip_count}")
    print(f"  Global minimum of partial sums: {min_partial:.8f} at n={min_partial_n}, K={min_partial_K}")

    # =========================================================================
    # 3. PHASE FORMULA ANALYSIS
    # =========================================================================
    print("\n--- Phase Formula Analysis ---")
    print("  lambda_n = 2 * sum_gamma [1 - cos(n * phi_gamma)]")
    print("  phi_gamma = pi - 2*arctan(2*gamma) in (0, pi)")

    zeros = load_zeros()
    phis = phase_from_gamma(zeros)

    print(f"\n  First 10 phase values phi_gamma = pi - 2*arctan(2*gamma):")
    print(f"  {'gamma':>10}  {'phi':>12}  {'1/gamma':>10}")
    for i in range(10):
        g = zeros[i]
        phi = phis[i]
        print(f"  {g:>10.4f}  {phi:>12.8f}  {1/g:>10.8f}")

    # Check: for small n, each term 2[1-cos(n*phi)] >= 0
    # When does the first term (gamma_1) start to decrease the sum?
    # It decreases when the term 2[1-cos(n*phi_1)] < 2[1-cos((n-1)*phi_1)]
    # i.e., cos(n*phi_1) > cos((n-1)*phi_1)
    # i.e., n*phi_1 > pi (first term decreases when n*phi_1 crosses pi)
    n_crossover_1 = np.pi / phis[0]
    print(f"\n  First zero crossover: n* = pi/phi_1 = pi/{phis[0]:.4f} = {n_crossover_1:.2f}")
    print(f"  => For n > {int(n_crossover_1)}, first zero's contribution starts oscillating down")

    # But is the TOTAL still increasing?
    # Compute lambda_n and lambda_{n+1} - lambda_n directly from 50K phases
    print(f"\n  Direct phase formula: lambda_n vs stored array:")
    print(f"  {'n':>6}  {'direct':>12}  {'stored':>12}  {'match':>6}")

    for n in [1, 2, 5, 10, 20, 44, 45, 46, 50, 100]:
        if n <= n_max:
            direct = 2.0 * np.sum(1.0 - np.cos(n * phis))
            stored = li[n-1]
            match = abs(direct - stored) < 1e-4
            print(f"  {n:>6}  {direct:>12.6f}  {stored:>12.6f}  {'YES' if match else 'NO':>6}")

    # Compute delta_n = lambda_{n+1} - lambda_n via phase formula
    print(f"\n  Delta = lambda_{{n+1}} - lambda_n via phase formula (n=40..55):")
    print(f"  (n_crossover = {n_crossover_1:.1f}, testing around there)")
    print(f"  {'n':>6}  {'delta':>12}  {'>0?':>5}")
    for n in range(40, 56):
        ln = 2.0 * np.sum(1.0 - np.cos(n * phis))
        ln1 = 2.0 * np.sum(1.0 - np.cos((n+1) * phis))
        delta = ln1 - ln
        print(f"  {n:>6}  {delta:>12.6f}  {'YES' if delta > 0 else 'NO':>5}")

    # =========================================================================
    # 4. FORMULA: lambda_{n+1} - lambda_n = 4 * sum_gamma sin((n+1/2)*phi)*sin(phi/2)
    # =========================================================================
    print("\n--- Analytical formula for delta_n ---")
    print("  delta_n = lambda_{n+1} - lambda_n")
    print("         = 4 * sum_gamma sin((n+1/2)*phi_gamma) * sin(phi_gamma/2)")
    print("  (Using: cos(A) - cos(B) = -2*sin((A+B)/2)*sin((A-B)/2))")
    print("  sin(phi_gamma/2) > 0 for all gamma (since 0 < phi < pi)")
    print("  But sin((n+1/2)*phi_gamma) can be negative for some gamma!")
    print()

    # Check: is delta_n > 0 for all n up to 50000?
    print("  Checking delta_n sign for all n=1..49999 via stored array...")
    neg_indices = np.where(deltas < 0)[0]
    if len(neg_indices) == 0:
        print("  RESULT: delta_n > 0 for ALL n=1..49999. Lambda is monotone increasing!")
        print(f"  Minimum delta: {min_delta:.8f} at n={min_delta_n}")
    else:
        print(f"  VIOLATIONS: delta_n < 0 at n = {neg_indices[:10] + 1} ...")

    # Histogram of deltas
    print(f"\n  Delta distribution:")
    print(f"  Min:    {np.min(deltas):.6f}")
    print(f"  Max:    {np.max(deltas):.6f}")
    print(f"  Mean:   {np.mean(deltas):.6f}")
    print(f"  Median: {np.median(deltas):.6f}")

    # =========================================================================
    # 5. MATHEMATICAL IMPLICATIONS
    # =========================================================================
    print("\n" + "=" * 70)
    print("MATHEMATICAL IMPLICATIONS")
    print("=" * 70)
    print("""
  MONOTONICITY OF LAMBDA_N:
  If lambda_n is strictly increasing for all n >= 1, then:
  - lambda_1 = p_1 = 1 + gamma_E/2 - log(2*sqrt(pi)) > 0 (unconditional!)
  - lambda_n >= lambda_1 > 0 for all n >= 1
  - This gives Li criterion: lambda_n > 0 for ALL n (under RH)

  IMPORTANT: Monotonicity holds NUMERICALLY for n=1..49999.
  But is it PROVABLE?

  Key formula: delta_n = 4 * sum_gamma sin((n+1/2)*phi_gamma) * sin(phi_gamma/2)

  Since sin(phi_gamma/2) > 0, the sign of delta_n depends on the sum
  of sin((n+1/2)*phi_gamma) weighted by sin(phi_gamma/2).

  This is a TRIGONOMETRIC SUM over the zeros. Proving it's positive
  for all n would be equivalent to proving RH (by Li's criterion).

  PARTIAL SUM STABILITY:
  If S_n(K) >= 0 for all K <= n, the Li coefficients are "monotone
  in the power sum truncation" — a structural property of the zeros.
  """)

    # Save results
    result = {
        "monotone_increasing": bool(neg_delta_count == 0),
        "min_delta": float(min_delta),
        "min_delta_n": int(min_delta_n),
        "neg_delta_count": int(neg_delta_count),
        "partial_sum_dip_count": dip_count,
        "global_min_partial": float(min_partial),
        "n_crossover_1": float(n_crossover_1),
    }
    out_path = os.path.join(RESULTS_DIR, "monotone_analysis.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"  Results saved to {out_path}")

if __name__ == "__main__":
    main()
