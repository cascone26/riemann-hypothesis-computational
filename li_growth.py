"""
Li coefficient growth rate analysis.

Known asymptotic (Keiper 1992, Li 1997):
  lambda_n ~ (n/2) * [log(n) + gamma - 1 - log(2*pi)] + O(log(n))

where gamma = Euler-Mascheroni constant.

If RH is true, all lambda_n > 0.
If RH is false, lambda_n eventually oscillates wildly.

We compare our computed coefficients against the asymptotic formula
and look for any anomalous behavior.
"""

import mpmath
import math
import json
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")


def load_zeros(count=None):
    path = os.path.join(DATA_DIR, "zeros_100k.txt")
    zeros = []
    with open(path, "r") as f:
        for i, line in enumerate(f):
            if count and i >= count:
                break
            line = line.strip()
            if line:
                zeros.append(float(line))
    return zeros


def compute_li_coefficients(zeros, n_max=200):
    """Compute Li coefficients using high-precision arithmetic."""
    mpmath.mp.dps = 50
    n_zeros = min(len(zeros), 20000)  # use 20K zeros for better convergence

    print(f"Computing Li coefficients lambda_1 to lambda_{n_max}")
    print(f"Using {n_zeros} zeros at 50 decimal places")
    print("-" * 60)

    results = []
    gamma = float(mpmath.euler)

    for n in range(1, n_max + 1):
        lambda_n = mpmath.mpf(0)
        for k in range(n_zeros):
            t = mpmath.mpf(zeros[k])
            rho = mpmath.mpc(0.5, t)
            term = 1 - (1 - 1 / rho) ** n
            rho_conj = mpmath.mpc(0.5, -t)
            term_conj = 1 - (1 - 1 / rho_conj) ** n
            lambda_n += term + term_conj

        lambda_n_float = float(lambda_n.real)

        # Keiper-Li asymptotic prediction
        if n >= 2:
            asymptotic = (n / 2) * (math.log(n) + gamma - 1 - math.log(2 * math.pi))
        else:
            asymptotic = 0.023  # known approximate value for n=1

        ratio = lambda_n_float / asymptotic if asymptotic != 0 else float('inf')

        results.append({
            "n": n,
            "lambda_n": lambda_n_float,
            "asymptotic": asymptotic,
            "ratio": ratio,
            "difference": lambda_n_float - asymptotic,
        })

        if n <= 20 or n % 20 == 0:
            print(f"  lambda_{n:>3d} = {lambda_n_float:>12.6f}  "
                  f"asymp = {asymptotic:>12.6f}  "
                  f"ratio = {ratio:.6f}")

    return results


def analyze_growth(results):
    """Analyze growth rate and look for anomalies."""
    print("\n" + "=" * 60)
    print("GROWTH ANALYSIS")
    print("=" * 60)

    # Check all positive
    all_positive = all(r["lambda_n"] > 0 for r in results)
    print(f"All lambda_n positive: {all_positive}")

    # Check ratio convergence (should approach 1 as n -> inf)
    last_10_ratios = [r["ratio"] for r in results[-10:]]
    mean_ratio = sum(last_10_ratios) / len(last_10_ratios)
    print(f"Mean ratio (last 10): {mean_ratio:.6f} (should approach 1.0)")

    # Check for any unexpected jumps
    for i in range(1, len(results)):
        prev = results[i - 1]["lambda_n"]
        curr = results[i]["lambda_n"]
        if prev > 0:
            growth = curr / prev
            if growth > 3.0 or growth < 0.5:
                print(f"  ANOMALY at n={results[i]['n']}: growth ratio {growth:.4f}")

    # Second differences (should be smooth for RH-consistent behavior)
    if len(results) >= 3:
        second_diffs = []
        for i in range(2, len(results)):
            d2 = results[i]["lambda_n"] - 2 * results[i-1]["lambda_n"] + results[i-2]["lambda_n"]
            second_diffs.append(d2)

        mean_d2 = sum(second_diffs) / len(second_diffs)
        max_d2 = max(abs(d) for d in second_diffs)
        print(f"Second differences: mean={mean_d2:.6f}, max_abs={max_d2:.6f}")
        print("  (Smooth second differences = consistent with RH)")

    return {"all_positive": all_positive, "asymptotic_ratio_convergence": mean_ratio}


def main():
    print("=" * 60)
    print("LI COEFFICIENT GROWTH RATE ANALYSIS")
    print("=" * 60)
    print()

    zeros = load_zeros(count=20000)
    print(f"Loaded {len(zeros)} high-precision zeros\n")

    results = compute_li_coefficients(zeros, n_max=200)
    analysis = analyze_growth(results)

    # Save
    outfile = os.path.join(RESULTS_DIR, "li_coefficients_200.json")
    with open(outfile, "w") as f:
        json.dump({"coefficients": results, "analysis": analysis}, f, indent=2)
    print(f"\nSaved to {outfile}")


if __name__ == "__main__":
    main()
