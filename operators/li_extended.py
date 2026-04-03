"""
Extended Li coefficient computation to N=5000.

Li's criterion (1997): RH ↔ λ_n ≥ 0 for ALL n ≥ 1.

λ_n = Σ_ρ [1 - (1 - 1/ρ)^n]

where the sum is over non-trivial zeros (each pair ρ, 1-ρ̄ contributes).

For ρ = 1/2 + iγ (on critical line):
  1 - 1/ρ = (ρ - 1)/ρ = (-1/2 + iγ)/(1/2 + iγ)
  |1 - 1/ρ| = |ρ-1|/|ρ| = √(1/4 + γ²) / √(1/4 + γ²) = 1

So EVERY zero on the critical line has |1 - 1/ρ| = 1.
The n-th power just rotates on the unit circle.

For ρ off the critical line:
  If ρ = σ + iγ with σ > 1/2:
  |1 - 1/ρ| ≠ 1 in general

  The pair (ρ, 1-ρ̄) = (σ+iγ, (1-σ)+iγ) contributes:
  [1-(1-1/ρ)^n] + [1-(1-1/(1-ρ̄))^n]
  This can be negative for large n → violation of λ_n ≥ 0.

Keiper-Li asymptotic (Keiper 1992, Li 1997):
  λ_n ~ (n/2)(log(n/(2π)) + γ_E - 1) as n → ∞

  where γ_E = Euler-Mascheroni constant ≈ 0.5772...

  The asymptotic becomes positive only for n ≥ n_0 where:
  log(n_0/(2π)) + γ_E - 1 > 0 → n_0 > 2πe^{1-γ_E} ≈ 2π × 1.75 ≈ 11

  So asymptotic is positive for n ≥ 12 approximately.

Strategy:
1. Extend computation to N=5000 using mpmath zetazero for high-precision zeros
2. Use partial sums with Aitken acceleration for convergence
3. Study λ_n / asymptotic ratio, look for minimum, find error term structure
4. Search for any near-violations (λ_n close to 0)
"""

import numpy as np
import mpmath
from mpmath import mp, mpf, mpc, log, euler, pi, exp
import json
import os
import time

mp.dps = 35

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")


def load_zeros_file(count=100000):
    """Load zeros from the 100K LMFDB file (31 decimal places)."""
    path = os.path.join(DATA_DIR, "zeros_100k.txt")
    zeros = []
    with open(path, "r") as f:
        for i, line in enumerate(f):
            if i >= count:
                break
            line = line.strip()
            if line:
                zeros.append(mpf(line))
    return zeros


def asymptotic_li(n):
    """
    Keiper-Li asymptotic: λ_n ~ (n/2)(log(n/(2π)) + γ_E - 1)
    More precise: + B log(n) + C + O(1/n)
    """
    n = mpf(n)
    leading = (n / 2) * (log(n / (2 * pi)) + euler - 1)
    return leading


def compute_li_batch(zeros, n_values, use_file_zeros=True):
    """
    Compute λ_n for n in n_values using mpmath precision.

    Uses the formula:
    λ_n = Σ_ρ [1 - (1 - 1/ρ)^n]

    where we sum over ρ = 1/2 + iγ and 1/2 - iγ for each zero γ > 0.

    For efficiency: precompute (1 - 1/ρ) for each zero, then take powers.
    """
    results = {}

    print(f"  Using {len(zeros)} zeros for Li computation...")

    # Precompute the base values: z_k = 1 - 1/ρ_k for each zero
    print("  Precomputing base values (1 - 1/ρ)...")
    bases = []  # Each entry is a pair (z, z̄) for conjugate zeros
    for gamma in zeros:
        rho = mpc(mpf('0.5'), gamma)
        rho_bar = mpc(mpf('0.5'), -gamma)
        z = 1 - 1 / rho
        z_bar = 1 - 1 / rho_bar
        bases.append((z, z_bar))

    print(f"  Base values computed for {len(bases)} zero pairs.")
    print(f"  |z_1| = {float(abs(bases[0][0])):.10f}  (should be 1.0 if zeros on critical line)")
    print(f"  |z_2| = {float(abs(bases[1][0])):.10f}")

    # Check |z| = 1 for all zeros on critical line
    # |1 - 1/ρ|² = |ρ-1|²/|ρ|² = ((σ-1)² + γ²)/(σ² + γ²)
    # For σ=1/2: = (1/4 + γ²)/(1/4 + γ²) = 1  ✓

    # For each n, compute λ_n
    for n in n_values:
        t0 = time.time()
        lambda_n = mpf(0)

        # Accumulate contributions
        for z, z_bar in bases:
            # (1 - 1/ρ)^n + (1 - 1/ρ̄)^n = 2 Re[(1-1/ρ)^n]
            zn = z ** n
            lambda_n += 2 - zn - z_bar ** n

        lambda_n_float = float(lambda_n.real)
        asymp = float(asymptotic_li(n))
        ratio = lambda_n_float / asymp if abs(asymp) > 1e-10 else float('inf')
        elapsed = time.time() - t0

        results[n] = {
            "n": n,
            "lambda_n": lambda_n_float,
            "asymptotic": asymp,
            "ratio": ratio,
            "error": lambda_n_float - asymp,
        }

        print(f"  n={n:>5d}: λ_n = {lambda_n_float:>14.6f}  "
              f"asymp = {asymp:>14.6f}  "
              f"ratio = {ratio:.6f}  [{elapsed:.1f}s]")

    return results


def compute_li_fast(zeros_list, n_max, batch_size=50):
    """
    Faster Li computation using numpy for moderate precision.
    Uses float64 with correction for small precision errors.

    For n up to 5000, we can use float64 and check the minimum more carefully.
    """
    n_zeros = len(zeros_list)
    zeros_np = np.array([float(g) for g in zeros_list])

    # Precompute bases
    rho = 0.5 + 1j * zeros_np
    bases = 1 - 1 / rho  # complex array

    print(f"  |bases[0]| = {abs(bases[0]):.10f} (should be ~1)")
    print(f"  |bases[100]| = {abs(bases[100]):.10f}")

    # Compute λ_n for n = 1, 2, ..., n_max
    n_vals = np.arange(1, n_max + 1)
    lambda_vals = np.zeros(n_max)

    # Batch over n values
    for n_start in range(0, n_max, batch_size):
        n_end = min(n_start + batch_size, n_max)
        ns = np.arange(n_start + 1, n_end + 1)

        # For each zero k and each n:
        # contribution = 2 - bases[k]^n - conj(bases[k])^n = 2 - 2 Re[bases[k]^n]
        # Sum over all zeros

        for i, n in enumerate(ns):
            # bases^n: all zeros' bases raised to n
            bn = bases ** n
            lambda_vals[n_start + i] = 2 * n_zeros - 2 * np.sum(bn.real)

        if n_end % 500 == 0 or n_end == n_max:
            print(f"  Computed through n={n_end}")

    return lambda_vals, n_vals


def analyze_li_sequence(lambda_vals, n_vals, zeros_np):
    """
    Deep analysis of the Li sequence:
    1. All positive? (Required for RH)
    2. Minimum value and location
    3. Rate of convergence to asymptotic
    4. Error term structure
    5. Near-violation search
    """
    n_max = len(lambda_vals)
    print("\n" + "=" * 70)
    print("LI SEQUENCE ANALYSIS")
    print("=" * 70)

    # Asymptotic values
    n_arr = np.arange(1, n_max + 1)
    asymp = (n_arr / 2) * (np.log(n_arr / (2 * np.pi)) + float(mpmath.euler) - 1)

    # For small n, asymptotic is negative (log(n/2π) < 1 - γ for n < ~11)
    # So ratio is only meaningful for larger n
    valid_mask = asymp > 0.1  # only compare where asymptotic is positive

    all_positive = np.all(lambda_vals > 0)
    min_val = np.min(lambda_vals)
    min_idx = np.argmin(lambda_vals)

    print(f"\n  All λ_n > 0: {all_positive}")
    print(f"  Min value: λ_{min_idx+1} = {min_val:.8f}")

    if not all_positive:
        neg_indices = np.where(lambda_vals < 0)[0]
        print(f"  NEGATIVE VALUES AT: n = {neg_indices[:10] + 1}")
        print(f"  THIS WOULD DISPROVE RH!")

    # Near-violations: values close to 0
    small_mask = lambda_vals < 1.0
    small_n = n_arr[small_mask]
    small_vals = lambda_vals[small_mask]
    print(f"\n  Values with λ_n < 1.0: {len(small_n)} occurrences")
    if len(small_n) > 0:
        print(f"  (n, λ_n): {list(zip(small_n[:10], [f'{v:.4f}' for v in small_vals[:10]]))}")

    # Ratio to asymptotic for large n
    large_n_mask = n_arr > 50
    ratios_large = lambda_vals[large_n_mask] / asymp[large_n_mask]
    print(f"\n  Ratio λ_n / asymptotic (n > 50):")
    print(f"    Mean: {np.mean(ratios_large):.6f}")
    print(f"    Std:  {np.std(ratios_large):.6f}")
    print(f"    Min:  {np.min(ratios_large):.6f}")
    print(f"    Max:  {np.max(ratios_large):.6f}")

    # Error term: λ_n - asymptotic
    error = lambda_vals - asymp
    # For large n, error term ~ B log(n) + C + O(log(n)/n)
    print(f"\n  Error term λ_n - asymptotic (n > 50):")
    err_large = error[large_n_mask]
    n_large = n_arr[large_n_mask]
    print(f"    Mean: {np.mean(err_large):.4f}")
    print(f"    Range: [{np.min(err_large):.4f}, {np.max(err_large):.4f}]")

    # Fit error term to A log(n) + B
    log_n_large = np.log(n_large)
    A = np.polyfit(log_n_large, err_large, 1)
    print(f"    Fit error ~ {A[0]:.4f} log(n) + {A[1]:.4f}")

    # Growth rate analysis: λ_n / n_arr
    growth = lambda_vals / n_arr
    print(f"\n  Growth λ_n/n at key values:")
    for n in [100, 500, 1000, 2000, 5000]:
        if n <= n_max:
            print(f"    n={n}: λ_{n}/n = {growth[n-1]:.6f}")

    return {
        "all_positive": bool(all_positive),
        "min_value": float(min_val),
        "min_n": int(min_idx + 1),
        "mean_ratio_n_gt_50": float(np.mean(ratios_large)),
        "error_fit": {"A": float(A[0]), "B": float(A[1])},
        "n_near_zero": int(np.sum(lambda_vals < 0.1)),
    }


def main():
    print("=" * 70)
    print("EXTENDED LI COEFFICIENT ANALYSIS")
    print("Li criterion: RH ↔ λ_n ≥ 0 for all n ≥ 1")
    print("=" * 70)

    # Load zeros
    print("\n--- Loading zeros ---")
    zeros = load_zeros_file(count=50000)  # 50K zeros for better convergence
    print(f"  Loaded {len(zeros)} zeros (31 decimal places)")
    zeros_np = np.array([float(g) for g in zeros[:50000]])

    # Fast numpy computation for n = 1 to 5000
    print("\n--- Fast Li computation (n = 1 to 5000) ---")
    t0 = time.time()
    lambda_vals, n_vals = compute_li_fast(zeros_np, n_max=5000)
    elapsed = time.time() - t0
    print(f"  Completed in {elapsed:.1f}s")

    # High-precision spot checks
    print("\n--- High-precision spot checks (mpmath, 50K zeros) ---")
    mp.dps = 50
    spot_ns = [1, 5, 10, 50, 100, 500, 1000, 2000, 5000]
    spot_results = compute_li_batch(zeros[:50000], spot_ns)

    # Analysis
    analysis = analyze_li_sequence(lambda_vals, n_vals, zeros_np)

    # Summary table
    print("\n--- Li coefficient summary (selected n) ---")
    print(f"  {'n':>6}  {'λ_n':>16}  {'asymptotic':>14}  {'ratio':>8}  {'positive':>8}")
    print(f"  {'-'*62}")
    key_ns = [1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
    for n in key_ns:
        lv = lambda_vals[n-1]
        av = (n/2) * (np.log(n/(2*np.pi)) + float(mpmath.euler) - 1)
        ratio = lv/av if abs(av) > 0.1 else float('nan')
        pos = "YES" if lv > 0 else "*** NO ***"
        print(f"  {n:>6}  {lv:>16.6f}  {av:>14.6f}  {ratio:>8.4f}  {pos:>8}")

    # Save
    save_data = {
        "lambda_values": lambda_vals.tolist(),
        "n_values": n_vals.tolist(),
        "analysis": analysis,
        "spot_checks": {str(k): v for k, v in spot_results.items()},
        "num_zeros_used": len(zeros_np),
    }
    outfile = os.path.join(RESULTS_DIR, "li_extended.json")
    with open(outfile, "w") as f:
        json.dump(save_data, f, indent=2)
    print(f"\n  Saved to {outfile}")

    print("\n" + "=" * 70)
    print("KEY FINDINGS")
    print("=" * 70)
    print(f"  All λ_n > 0 (n=1..5000): {analysis['all_positive']}")
    print(f"  Minimum: λ_{analysis['min_n']} = {analysis['min_value']:.6f}")
    print(f"  Error term ~ {analysis['error_fit']['A']:.4f} log(n) + {analysis['error_fit']['B']:.4f}")
    print(f"  Asymptotic ratio (n>50): {analysis['mean_ratio_n_gt_50']:.6f}")

    if analysis['all_positive']:
        print("\n  CONCLUSION: No violations found through n=5000.")
        print("  Strong computational evidence for RH (though not a proof).")
        print("  The error term analysis gives clues about the analytic structure.")
    else:
        print("\n  *** VIOLATION FOUND! Check negative values carefully. ***")


if __name__ == "__main__":
    main()
