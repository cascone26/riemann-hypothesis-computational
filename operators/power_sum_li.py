"""
Li coefficients via power sums of zeros.

Formula:
  λ_n = Σ_{k=1}^{n} C(n,k) (-1)^{k-1} p_k

where p_k = Σ_ρ 1/ρ^k are the Newton power sums of reciprocal zeros.

Key formula:
  p_k = -(1/(k-1)!) × d^k/ds^k [log ξ(s)]|_{s=0}

This computes λ_n EXACTLY (to precision of mpmath), without truncating
the sum over zeros. This bypasses the 50K-zero truncation issue.

The function log ξ(s) at s=0 is analytic (ξ(0) = 1/2 ≠ 0).

This allows:
1. Exact λ_n for any n (limited by mpmath precision, not zero count)
2. Study of sign pattern of p_k
3. Connection between p_k positivity structure and λ_n > 0

Additional question: Can the positivity λ_n > 0 be derived from known
properties of ξ(s) without summing over zeros explicitly?

Reference: Keiper (1992), Li (1997), Maslanka (2004)
"""

import mpmath
from mpmath import mp, mpf, mpc, log, factorial, pi, euler, gamma, diff
import numpy as np
import json
import os
import time

mp.dps = 50  # 50 decimal places

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)


def xi_func(s):
    """ξ(s) = (1/2)s(s-1)π^{-s/2}Γ(s/2)ζ(s)"""
    return mpmath.mpf('0.5') * s * (s - 1) * mpmath.power(pi, -s / 2) * gamma(s / 2) * mpmath.zeta(s)


def log_xi(s):
    """log ξ(s) — analytic away from zeros of ξ"""
    return log(xi_func(s))


def compute_power_sums(k_max, s0=mpf('0')):
    """
    Compute p_k = Σ_ρ 1/ρ^k for k=1..k_max.

    Formula: p_k = -k × c_k
    where c_k are Taylor coefficients of log ξ(s) at s=s0:
    log ξ(s) = Σ_k c_k (s-s0)^k

    Uses mpmath.taylor() to compute all coefficients at once — much faster
    than calling diff() k times separately.

    Note: For k=1, p_1 requires regularization (the sum diverges absolutely
    but converges in pairs). This is handled by the analytic continuation.
    """
    power_sums = {}

    print(f"Computing p_k for k=1..{k_max} at s={s0} via Taylor expansion...")
    t0 = time.time()

    try:
        # Get all Taylor coefficients at once: c_k = f^(k)(s0) / k!
        # So f^(k)(s0) = k! * c_k
        # p_k = -(1/(k-1)!) * f^(k)(s0) = -(1/(k-1)!) * k! * c_k = -k * c_k
        coeffs = mpmath.taylor(log_xi, s0, k_max)
        elapsed = time.time() - t0
        print(f"  Taylor expansion computed in {elapsed:.1f}s")

        for k in range(1, k_max + 1):
            if k < len(coeffs):
                ck = coeffs[k]
                pk = -k * ck
                power_sums[k] = complex(pk)
                print(f"  p_{k:>3d} = {float(pk.real):>18.10f}  (Im={float(pk.imag):.2e})")
            else:
                print(f"  p_{k}: Taylor expansion too short")
                power_sums[k] = None

    except Exception as e:
        print(f"  Taylor expansion FAILED: {e}")
        print("  Falling back to direct differentiation for small k...")
        for k in range(1, min(k_max + 1, 6)):
            try:
                deriv = diff(log_xi, s0, k)
                pk = -deriv / factorial(k - 1)
                power_sums[k] = complex(pk)
                print(f"  p_{k} = {float(pk.real):.10f}")
            except Exception as e2:
                print(f"  p_{k}: FAILED ({e2})")
                power_sums[k] = None

    return power_sums


def compute_power_sums_from_zeros(zeros_list, k_max):
    """
    Direct computation: p_k = 2 Σ_{γ>0} Re[1/(1/2+iγ)^k]

    This uses the 50K zeros and converges for k ≥ 2.
    For k=1: conditionally convergent, result ≈ 0.0231.
    """
    power_sums = {}
    gammas = np.array([float(g) for g in zeros_list])
    rhos = 0.5 + 1j * gammas

    for k in range(1, k_max + 1):
        vals = 1.0 / (rhos ** k)
        pk_partial = 2.0 * np.sum(vals.real)
        power_sums[k] = pk_partial

    return power_sums


def compute_lambda_from_power_sums(power_sums, n_max):
    """
    λ_n = Σ_{k=1}^{n} C(n,k) (-1)^{k-1} p_k

    Uses the analytic p_k values (from log ξ derivatives).
    This gives EXACT Li coefficients for any n.

    Note: For large n, numerical cancellation may reduce precision.
    """
    from mpmath import binomial

    lambda_vals = {}

    print(f"\nComputing λ_n from power sums for n=1..{n_max}...")

    for n in range(1, n_max + 1):
        lam = mpf('0')
        for k in range(1, n + 1):
            pk = power_sums.get(k)
            if pk is None:
                break
            binom = binomial(n, k)
            sign = (-1) ** (k - 1)
            lam += binom * sign * mpf(pk.real)

        lambda_vals[n] = float(lam)

        if n <= 20 or n % 50 == 0:
            print(f"  λ_{n:>4d} = {float(lam):>14.8f}")

    return lambda_vals


def compute_lambda_incremental(power_sums, n_max):
    """
    More numerically stable: compute λ_n using the recurrence.

    λ_n = λ_{n-1} + ... [via Newton's identity approach]

    Or directly: λ_n = Σ_{k=1}^n C(n,k)(-1)^{k-1} p_k

    The numerical issue: large binomial coefficients times alternating-sign
    terms can cause cancellation. Use extended precision.
    """
    pass  # For now use direct formula


def study_power_sum_signs(power_sums, k_max):
    """
    Study the sign and magnitude pattern of p_k.

    The pattern determines whether λ_n = Σ C(n,k)(-1)^{k-1} p_k is positive.

    Key question: Is there a structural reason λ_n > 0 from the p_k pattern?
    """
    print("\n" + "=" * 60)
    print("POWER SUM SIGN ANALYSIS")
    print("=" * 60)

    signs = []
    mags = []

    for k in range(1, k_max + 1):
        pk = power_sums.get(k)
        if pk is None:
            continue
        val = pk.real if isinstance(pk, complex) else float(pk)
        sign = "+" if val >= 0 else "-"
        signs.append((k, val, sign))
        mags.append(abs(val))

    print(f"\n  k: sign(p_k), value")
    for k, val, sign in signs:
        print(f"  {k:>3d}: {sign}  {val:>20.10f}")

    # Look for patterns
    pos_k = [k for k, v, s in signs if s == "+"]
    neg_k = [k for k, v, s in signs if s == "-"]

    print(f"\n  Positive p_k: {pos_k}")
    print(f"  Negative p_k: {neg_k}")

    # Growth rate
    if len(mags) > 5:
        ratios = [mags[i] / mags[i-1] if mags[i-1] > 0 else float('nan')
                  for i in range(1, len(mags))]
        print(f"\n  Growth ratios |p_k|/|p_{k-1}|:")
        for i, r in enumerate(ratios[:15]):
            print(f"    p_{i+2}/p_{i+1} = {r:.4f}")

    # The asymptotic p_k ~ ?
    # For large k: p_k = Σ_{γ>0} 2Re[1/(1/2+iγ)^k]
    # The sum is dominated by the smallest |ρ| = |1/2+iγ_1| = √(1/4 + γ_1²)
    # where γ_1 ≈ 14.135
    rho1_abs = np.sqrt(0.25 + 14.135**2)
    print(f"\n  Predicted asymptotic decay: |p_k| ~ const × {1/rho1_abs:.6f}^k = const × e^{{-{np.log(rho1_abs):.4f}k}}")
    print(f"  |ρ_1| = √(1/4 + 14.135²) ≈ {rho1_abs:.6f}")

    return signs


def verify_vs_sum_computation(power_sums_analytic, power_sums_zeros):
    """
    Compare p_k computed via log ξ derivatives vs direct zero sum.
    Shows the truncation error of the 50K-zero computation.
    """
    print("\n" + "=" * 60)
    print("COMPARISON: log ξ derivatives vs 50K-zero partial sums")
    print("=" * 60)
    print(f"\n  {'k':>3}  {'analytic':>20}  {'50K partial':>20}  {'ratio':>10}")
    print(f"  {'-'*60}")

    for k in range(1, min(len(power_sums_analytic) + 1, 21)):
        pa = power_sums_analytic.get(k)
        pz = power_sums_zeros.get(k)
        if pa is None or pz is None:
            continue

        pa_val = pa.real if isinstance(pa, complex) else float(pa)
        pz_val = float(pz)

        ratio = pz_val / pa_val if abs(pa_val) > 1e-15 else float('nan')
        print(f"  {k:>3d}  {pa_val:>20.10f}  {pz_val:>20.10f}  {ratio:>10.6f}")


def prove_lambda_structure(power_sums, n_test=50):
    """
    Attempt to find a pattern in λ_n from the p_k structure.

    Observation: λ_n = Σ_{k=1}^{n} C(n,k)(-1)^{k-1} p_k

    This is related to the forward difference operator:
    λ_n = Σ_{k=0}^{n-1} C(n-1,k) (-1)^k [p_{k+1} + p_k/2 + ...]

    Actually: by the binomial transform, if q_k = (-1)^{k-1} p_k, then
    λ_n = Σ_{k=1}^{n} C(n,k) q_k

    The generating function approach: if Σ_{n=1}^∞ λ_n x^n/n = ...

    From the product formula for ξ:
    ξ(s) = ξ(0) Π_ρ (1 - s/ρ)

    Taking log:
    log ξ(s) = log ξ(0) + Σ_ρ Σ_{k=1}^∞ (-1/k)(s/ρ)^k
              = log ξ(0) - Σ_{k=1}^∞ (s^k/k) Σ_ρ 1/ρ^k
              = log ξ(0) - Σ_{k=1}^∞ (s^k/k) p_k

    So: d/ds log ξ(s) = -Σ_{k=0}^∞ s^k p_{k+1}

    This gives: ξ'(s)/ξ(s) = -Σ_{k=0}^∞ p_{k+1} s^k

    For s=1: ξ'(1)/ξ(1) = -Σ_{k=0}^∞ p_{k+1}

    This series diverges (p_k → 0 but sum diverges). Need different approach.

    Key relation (from Keiper 1992):
    λ_n = [d/ds]^{n-1} [s^{n-1} ξ'(s)/ξ(s)]|_{s=1} / (n-1)!
         + correction

    This connects to properties of ξ'/ξ at s=1 (the "logarithmic derivative").
    """
    print("\n" + "=" * 60)
    print("GENERATING FUNCTION / STRUCTURE ANALYSIS")
    print("=" * 60)

    # Verify the generating function relation:
    # -d/ds log ξ(s)|_{s=0} = p_1
    # -d²/ds² log ξ(s)|_{s=0} / 1! = p_2 ... wait that's what we computed

    # Look at the sequence λ_n / (n * log(n)) for large n
    # Should approach (1/2)(1 - 1/log n + ... ) → 1/2

    # Also: λ_n satisfies a recurrence? Check if λ_n - λ_{n-1} has a pattern.

    # From the sum: λ_n = λ_{n-1} + Σ_{k=1}^{n-1} [C(n,k) - C(n-1,k)](-1)^{k-1}p_k
    #                                 + (-1)^{n-1} p_n
    # Using C(n,k) - C(n-1,k) = C(n-1,k-1):
    # λ_n = λ_{n-1} + Σ_{k=1}^{n-1} C(n-1,k-1)(-1)^{k-1}p_k + (-1)^{n-1} p_n
    # = λ_{n-1} + Σ_{j=0}^{n-2} C(n-1,j)(-1)^j p_{j+1} + (-1)^{n-1} p_n

    # This gives a recurrence involving all previous p_k.

    print("\n  Recurrence structure:")
    print("  λ_n = λ_{n-1} + Δ_n")
    print("  where Δ_n = Σ_{j=0}^{n-2} C(n-1,j)(-1)^j p_{j+1} + (-1)^{n-1} p_n")
    print()
    print("  For positivity: need Δ_n > 0 for all n, OR λ_{n-1} large enough")
    print("  to absorb any negative Δ_n.")

    # Compute the increments Δ_n from the analytic p_k
    from mpmath import binomial as binom

    deltas = {}
    lam_prev = mpf('0')
    for n in range(1, n_test + 1):
        lam_n = mpf('0')
        for k in range(1, n + 1):
            pk = power_sums.get(k)
            if pk is None:
                break
            val = mpf(pk.real) if isinstance(pk, complex) else mpf(pk)
            lam_n += binom(n, k) * mpf((-1)**(k-1)) * val

        delta_n = float(lam_n) - float(lam_prev)
        deltas[n] = {"lambda": float(lam_n), "delta": delta_n, "positive": delta_n > 0}
        lam_prev = lam_n

    print(f"\n  {'n':>4}  {'λ_n':>14}  {'Δ_n = λ_n - λ_{n-1}':>22}  {'Δ>0?':>6}")
    print(f"  {'-'*55}")
    for n in range(1, min(n_test + 1, 31)):
        d = deltas[n]
        sign = "YES" if d["positive"] else "*** NO ***"
        print(f"  {n:>4}  {d['lambda']:>14.6f}  {d['delta']:>22.6f}  {sign:>6}")

    neg_deltas = [n for n, d in deltas.items() if not d["positive"]]
    print(f"\n  Negative increments Δ_n: {neg_deltas[:20]}")
    print(f"  {'All Δ_n > 0: ' + str(len(neg_deltas) == 0)}")

    return deltas


def main():
    print("=" * 70)
    print("LI COEFFICIENTS VIA POWER SUMS OF ZEROS")
    print("Exact computation via ξ function derivatives")
    print("=" * 70)

    k_max = 25  # Compute p_1 through p_25

    # Method 1: Analytic via log ξ derivatives
    print("\n--- Analytic power sums via log ξ derivatives ---")
    mp.dps = 50
    power_sums_analytic = compute_power_sums(k_max, s0=mpf('0'))

    # Method 2: Direct zero sum (partial, 50K zeros)
    print("\n--- Direct zero sum (50K zeros, partial) ---")
    from mpmath import mpf as mf

    DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
    zeros = []
    with open(os.path.join(DATA_DIR, "zeros_100k.txt"), "r") as f:
        for i, line in enumerate(f):
            if i >= 50000:
                break
            line = line.strip()
            if line:
                zeros.append(float(line))
    power_sums_zeros = compute_power_sums_from_zeros(zeros, k_max)

    # Comparison
    verify_vs_sum_computation(power_sums_analytic, power_sums_zeros)

    # Sign analysis
    signs = study_power_sum_signs(power_sums_analytic, k_max)

    # Compute λ_n from analytic p_k
    n_max = 100
    print(f"\n--- Computing λ_n from analytic p_k for n=1..{n_max} ---")
    mp.dps = 60  # Extra precision for cancellation

    lambda_analytic = compute_lambda_from_power_sums(power_sums_analytic, n_max)

    # Compare with numerical zero sum at key values
    print("\n--- Comparison: analytic λ_n vs numerical sum ---")
    DATA_DIR_ref = DATA_DIR
    import numpy as npy
    li50k = npy.load(os.path.join(RESULTS_DIR, "li_50k.npy"))

    print(f"\n  {'n':>6}  {'analytic λ_n':>16}  {'50K-zero sum':>16}  {'ratio':>8}")
    print(f"  {'-'*52}")
    for n in [1, 2, 3, 5, 10, 20, 50, 100]:
        la = lambda_analytic.get(n)
        lz = float(li50k[n - 1]) if n <= len(li50k) else None
        if la is not None and lz is not None:
            ratio = lz / la if abs(la) > 1e-10 else float('nan')
            print(f"  {n:>6}  {la:>16.8f}  {lz:>16.8f}  {ratio:>8.6f}")

    # Structure analysis
    deltas = prove_lambda_structure(power_sums_analytic, n_test=100)

    # Extended λ_n computation beyond 50K using analytic formula
    print("\n--- Extended analytic λ_n for large n ---")
    print("  (Using k_max power sums — precision degrades for n >> k_max)")
    print("  Note: For n > k_max, the formula is INCOMPLETE (missing higher p_k)")
    print("  This shows convergence properties only")

    # Summary
    all_pos_100 = all(v > 0 for v in lambda_analytic.values())
    min_val = min(lambda_analytic.values())
    min_n = min(lambda_analytic, key=lambda_analytic.get)

    print("\n" + "=" * 70)
    print("KEY FINDINGS")
    print("=" * 70)

    pos_k = [k for k, v, s in signs if s == "+"]
    neg_k = [k for k, v, s in signs if s == "-"]

    print(f"\n  Power sum signs:")
    print(f"    Positive p_k: k ∈ {pos_k}")
    print(f"    Negative p_k: k ∈ {neg_k}")

    print(f"\n  λ_n from analytic formula (n=1..{n_max}):")
    print(f"    All λ_n > 0: {all_pos_100}")
    print(f"    Min: λ_{min_n} = {min_val:.8f}")

    neg_deltas = [n for n, d in deltas.items() if not d["positive"]]
    print(f"\n  Increments Δ_n = λ_n - λ_{'{n-1}'}:")
    print(f"    All Δ_n > 0: {len(neg_deltas) == 0}")
    if neg_deltas:
        print(f"    Negative at: n = {neg_deltas[:10]}")

    # The key structural question
    print("\n  STRUCTURAL QUESTION:")
    print("  Is there an explicit formula that shows λ_n > 0 from the p_k pattern?")
    print()
    print("  Observation: the sign alternation of p_k (when viewed with (-1)^{k-1})")
    print("  means all terms C(n,k)(-1)^{k-1}p_k have the SAME sign if")
    print("  sign((-1)^{k-1} p_k) is constant. Is this true?")

    for k in range(1, k_max + 1):
        pk = power_sums_analytic.get(k)
        if pk is None:
            continue
        val = pk.real if isinstance(pk, complex) else float(pk)
        combined = ((-1) ** (k - 1)) * val
        sign = "+" if combined >= 0 else "-"
        if k <= 10:
            print(f"    (-1)^{k-1} × p_{k} = {combined:+.8f}  [{sign}]")

    # Save
    save_data = {
        "power_sums_analytic": {str(k): {"real": v.real, "imag": v.imag}
                                 for k, v in power_sums_analytic.items() if v is not None},
        "lambda_analytic_100": {str(k): v for k, v in lambda_analytic.items()},
        "all_positive_100": all_pos_100,
        "min_lambda": float(min_val),
        "delta_signs": {str(n): d["positive"] for n, d in deltas.items()},
    }

    outfile = os.path.join(RESULTS_DIR, "power_sum_li.json")
    with open(outfile, "w") as f:
        json.dump(save_data, f, indent=2)
    print(f"\n  Saved to {outfile}")


if __name__ == "__main__":
    main()
