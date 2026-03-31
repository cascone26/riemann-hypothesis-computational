"""
Computational verification of Robin's inequality:
  sigma_1(n) < e^gamma * n * ln(ln(n))  for n > 5040

We check this for all n from 5041 to some large bound, identifying:
1. How close the inequality gets to being tight
2. The "worst offenders" (largest sigma(n) / robin_bound(n) ratio)
3. Whether there are any violations (there shouldn't be if RH is true)

5040 = 7! is the last integer where the inequality FAILS.
It's a "colossally abundant number" — has unusually large sigma(n)/n.
"""

import math
import json
import os

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")

def sigma1(n):
    """Sum of divisors of n."""
    if n <= 0:
        return 0
    total = 0
    for d in range(1, int(math.isqrt(n)) + 1):
        if n % d == 0:
            total += d
            if d != n // d:
                total += n // d
    return total


def robin_bound(n):
    """e^gamma * n * ln(ln(n))"""
    gamma = 0.5772156649015329
    if n <= 2:
        return float('inf')
    log_n = math.log(n)
    if log_n <= 0:
        return float('inf')
    log_log_n = math.log(log_n)
    return math.exp(gamma) * n * log_log_n


def verify_robin(n_max=100000):
    print(f"Verifying Robin's inequality for n = 5041 to {n_max}")
    print("=" * 60)

    violations = []
    worst_ratio = 0
    worst_n = 0
    ratios_near_1 = []

    # Check n = 5040 itself (should FAIL)
    s = sigma1(5040)
    b = robin_bound(5040)
    print(f"\nn = 5040 (should FAIL):")
    print(f"  sigma(5040) = {s}")
    print(f"  robin_bound = {b:.4f}")
    print(f"  ratio = {s/b:.8f}")
    print(f"  Fails: {s >= b}")

    print(f"\nChecking n = 5041 to {n_max}...")

    for n in range(5041, n_max + 1):
        s = sigma1(n)
        b = robin_bound(n)
        ratio = s / b

        if ratio >= 1.0:
            violations.append((n, s, b, ratio))

        if ratio > worst_ratio:
            worst_ratio = ratio
            worst_n = n

        if ratio > 0.95:
            ratios_near_1.append((n, ratio))

    print(f"\nResults:")
    print(f"  Violations: {len(violations)}")
    print(f"  Worst ratio: {worst_ratio:.8f} at n = {worst_n}")
    print(f"  sigma({worst_n}) = {sigma1(worst_n)}")
    print(f"  robin_bound({worst_n}) = {robin_bound(worst_n):.4f}")

    if violations:
        print(f"\n  VIOLATIONS FOUND (would disprove RH!):")
        for n, s, b, r in violations[:10]:
            print(f"    n={n}: sigma={s}, bound={b:.4f}, ratio={r:.8f}")
    else:
        print(f"\n  NO VIOLATIONS — consistent with RH")

    # Show the closest calls
    ratios_near_1.sort(key=lambda x: -x[1])
    print(f"\n  Closest to violating (ratio > 0.95):")
    for n, r in ratios_near_1[:20]:
        print(f"    n = {n:>8d}: ratio = {r:.8f}")

    # Highly composite numbers are the interesting ones
    # Check some known highly composite numbers
    hcn = [5040, 7560, 10080, 15120, 20160, 25200, 27720, 45360,
           50400, 55440, 83160, 110880, 166320, 221760, 277200,
           332640, 498960, 554400, 665280]

    print(f"\n  Highly composite numbers:")
    for n in hcn:
        if n > 5040:
            s = sigma1(n)
            b = robin_bound(n)
            print(f"    n = {n:>8d}: sigma = {s:>10d}, bound = {b:>12.2f}, ratio = {s/b:.8f}")

    results = {
        "n_max": n_max,
        "violations": len(violations),
        "worst_ratio": worst_ratio,
        "worst_n": worst_n,
        "n_checked": n_max - 5040,
    }
    with open(os.path.join(RESULTS_DIR, "robin_verification.json"), "w") as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    verify_robin(100000)
