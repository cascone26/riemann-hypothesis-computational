"""
Jensen Polynomial Hyperbolicity — Newton Recurrence Method

Instead of mpmath.diff (slow), use Newton-Girard recurrence:
  GORZ_a_n = (1/n) Σ_{k=1}^n (-1)^{k+1} P_{2k} GORZ_a_{n-k}

where P_{2k} = Σ_{γ>0} 1/γ^{2k} (power sums of reciprocals of zeta zeros squared).

GORZ notation: Ξ(z) = Σ_{n≥0} (-1)^n GORZ_a_n z^{2n}, GORZ_a_n > 0 for all n.
Substituting w = -z²: F(w) = Σ GORZ_a_n w^n = Ξ(√(-w)).
Euler product: F(w) = Ξ(0) × Π_{γ>0} (1 + w/γ²).
log F = log Ξ(0) + Σ_{k≥1} (-1)^{k+1} P_{2k} w^k / k.
Newton: n GORZ_a_n = Σ_{k=1}^n (-1)^{k+1} P_{2k} GORZ_a_{n-k}.

Griffin-Ono-Rolen-Zagier (2019): J_{d,n} hyperbolic for d≤8, all n.
Goal: verify d=9..30, n=0..200.
"""

import mpmath
import numpy as np
from math import comb
import json
import os
import time

mpmath.mp.dps = 80  # 80 decimal places

ZERO_FILE   = "/Users/scones/projects/riemann/data/zeros_odlyzko_100k"
RESULTS_DIR = "/Users/scones/projects/riemann/results"

# ── Load zeros ───────────────────────────────────────────────────────────────

def load_zeros(n_max=20000):
    """Load first n_max zeros from Odlyzko file."""
    gammas = []
    with open(ZERO_FILE) as f:
        for line in f:
            line = line.strip()
            if line:
                gammas.append(float(line))
                if len(gammas) >= n_max:
                    break
    return gammas

# ── Compute power sums ────────────────────────────────────────────────────────

def compute_power_sums(gammas, k_max, add_tail=True):
    """
    P_{2k} = Σ_{γ>0} 1/γ^{2k}  for k = 1..k_max.

    Uses mpmath for accuracy. Adds analytic tail estimate for k=1,2.
    For k≥3, first 10k zeros give ~machine precision.

    Returns list P[k] for k=0..k_max (P[0] = # zeros, P[k] = Σ 1/γ^{2k}).
    """
    print(f"Computing P_{{2k}} for k=1..{k_max} using {len(gammas)} zeros...")
    t0 = time.time()

    # Use mpmath for the full computation
    mp_gammas = [mpmath.mpf(str(g)) for g in gammas]

    P = [mpmath.mpf('0')] * (k_max + 1)
    P[0] = mpmath.mpf(len(gammas))

    for k in range(1, k_max + 1):
        # For large k, only small γ contribute (1/γ^{2k} underflows for large γ)
        # Cutoff: if γ^{2k} > 10^{k_max}, skip. γ_cutoff = 10^{k_max/(2k)}.
        # For k=1: cutoff ≈ 10^{k_max/2} ≈ 10^{112} (all zeros count)
        # For k=k_max: cutoff ≈ 10^{0.5} ≈ 3.16 (only first zero contributes)

        p2k = mpmath.mpf('0')
        for g in mp_gammas:
            p2k += mpmath.power(g, -2*k)
        P[k] = p2k

        if (k % 25 == 0) or k == k_max:
            elapsed = time.time() - t0
            print(f"  k={k}: P_{{2k}} = {float(P[k]):.6e}  ({elapsed:.1f}s)")

    # Tail correction for k=1 (only significant case):
    # Σ_{γ > γ_N} 1/γ^2 ≈ ∫_N^∞ (log(t)/(2π))/t^2 dt ≈ (log N + 1) / (2π N)
    # where N = γ_{last used zero} (index-based cutoff via N(T) ~ T log T / 2π)
    # Actually: approximate N(T) ≈ T/(2π) × log(T/(2πe))
    # Use simpler: tail ≈ (log γ_N + 1) / (2π × γ_N) for k=1
    if add_tail and len(gammas) >= 100:
        g_N = mpmath.mpf(str(gammas[-1]))
        tail_k1 = (mpmath.log(g_N) + 1) / (2 * mpmath.pi * g_N)
        old = P[1]
        P[1] += tail_k1
        print(f"  Tail correction for k=1: {float(tail_k1):.4e} (was {float(old):.6e}, now {float(P[1]):.6e})")

    print(f"Power sums done in {time.time()-t0:.1f}s")
    return P

# ── Newton recurrence for GORZ a_n ───────────────────────────────────────────

def compute_gorz_coefficients(P, n_max, Xi0=None):
    """
    Newton-Girard recurrence:
      GORZ_a_n = (1/n) Σ_{k=1}^n (-1)^{k+1} P_{2k} GORZ_a_{n-k}

    GORZ_a_0 = Ξ(0) = ξ(1/2).

    Returns list a[0..n_max] as mpmath numbers.
    """
    if Xi0 is None:
        # Ξ(0) = ξ(1/2) = (1/2)(1/2)(-1/2) π^{-1/4} Γ(1/4) ζ(1/2)
        s = mpmath.mpf('0.5')
        Xi0 = mpmath.re(mpmath.mpf('0.5') * s * (s-1) *
                        mpmath.power(mpmath.pi, -s/2) *
                        mpmath.gamma(s/2) * mpmath.zeta(s))
        print(f"Ξ(0) = {float(Xi0):.10f}")

    a = [mpmath.mpf('0')] * (n_max + 1)
    a[0] = mpmath.mpf(str(float(Xi0)))

    print(f"\nComputing GORZ a_n for n=0..{n_max} via Newton recurrence...")
    t0 = time.time()

    for n in range(1, n_max + 1):
        s = mpmath.mpf('0')
        for k in range(1, n + 1):
            sign = mpmath.mpf('1') if (k % 2 == 1) else mpmath.mpf('-1')
            s += sign * P[k] * a[n - k]
        a[n] = s / n

        if n <= 10 or n % 25 == 0 or n == n_max:
            print(f"  a[{n:3d}] = {float(a[n]):.6e}  ({time.time()-t0:.1f}s)")

    print(f"Coefficients done in {time.time()-t0:.1f}s")
    return a

# ── Jensen polynomial hyperbolicity ──────────────────────────────────────────

def jensen_poly_coeffs(a, d, n):
    """
    J_{d,n}(X) = Σ_{j=0}^d C(d,j) GORZ_a_{n+j} X^j

    Returns coefficients [c_0, c_1, ..., c_d] as floats.
    (c_j = C(d,j) × a[n+j])
    """
    return [float(comb(d, j)) * float(a[n + j]) for j in range(d + 1)]

def is_hyperbolic(coeffs, tol=1e-8):
    """
    Check if polynomial (given by [c_0, c_1, ..., c_d]) is hyperbolic (all real roots).
    Uses numpy.roots on the reversed (high-degree-first) form.
    """
    # numpy.roots expects [a_d, a_{d-1}, ..., a_0]
    p = np.array(coeffs[::-1], dtype=complex)
    if len(p) <= 1:
        return True, 0.0
    roots = np.roots(p)
    max_imag = max(abs(r.imag) for r in roots) if len(roots) > 0 else 0.0
    scale    = max(abs(r.real) for r in roots) + 1e-30 if len(roots) > 0 else 1.0
    return bool(max_imag / scale < tol), float(max_imag)

def run_hyperbolicity_table(a, d_max, n_max):
    """
    For each d from 1 to d_max, check J_{d,n} for n=0..min(n_max, len(a)-d-1).
    Returns: dict {d: {"n_tested": ..., "n_hyperbolic": ..., "first_fail": ...}}
    """
    results = {}
    print(f"\n{'='*70}")
    print(f"JENSEN POLYNOMIAL HYPERBOLICITY TABLE")
    print(f"d = 1..{d_max}, n = 0..{n_max}")
    print(f"{'='*70}")
    print(f"{'d':>4}  {'n_tested':>8}  {'n_hyp':>6}  {'first_fail':>10}  {'min_ratio':>12}")
    print(f"{'-'*50}")

    for d in range(1, d_max + 1):
        n_available = len(a) - d - 1
        n_limit = min(n_max, n_available)
        if n_limit < 0:
            print(f"{d:>4}  not enough coefficients")
            break

        n_tested = 0
        n_hyp    = 0
        first_fail = None
        min_ratio  = float('inf')  # min(max_imag / scale) — lower = more hyperbolic

        for n in range(n_limit + 1):
            coeffs = jensen_poly_coeffs(a, d, n)
            hyp, max_imag = is_hyperbolic(coeffs)
            scale = max(abs(c) for c in coeffs) + 1e-30
            ratio = max_imag / scale

            n_tested += 1
            if hyp:
                n_hyp += 1
            else:
                if first_fail is None:
                    first_fail = n
            if ratio < min_ratio:
                min_ratio = ratio

        status = "✓ ALL" if n_hyp == n_tested else f"✗ FAIL@n={first_fail}"
        print(f"{d:>4}  {n_tested:>8}  {n_hyp:>6}  {str(first_fail):>10}  {min_ratio:>12.4e}  {status}")

        results[d] = {
            "d": d,
            "n_tested": n_tested,
            "n_hyperbolic": n_hyp,
            "first_fail": first_fail,
            "min_ratio": min_ratio,
            "all_hyperbolic": n_hyp == n_tested
        }

    return results

# ── Discriminant analysis ─────────────────────────────────────────────────────

def analyze_discriminant_growth(a, d, n_vals):
    """
    For fixed d, compute discriminant of J_{d,n} for various n.
    Discriminant > 0 ↔ all real roots (for odd degree: discriminant has specific sign).
    Tracks the MINIMUM imaginary part ratio across tested n.
    """
    print(f"\nDiscriminant analysis for d={d}:")
    for n in n_vals:
        if n + d >= len(a):
            break
        coeffs = jensen_poly_coeffs(a, d, n)
        hyp, max_imag = is_hyperbolic(coeffs, tol=1e-12)
        roots = np.roots(np.array(coeffs[::-1], dtype=complex))
        root_info = [(float(r.real), float(r.imag)) for r in roots]
        root_info.sort(key=lambda x: x[0])
        print(f"  n={n:4d}: hyp={'✓'if hyp else '✗'}  max_imag={max_imag:.3e}  roots: {root_info[:3]}...")

# ── Log-convexity check (needed for d=2 hyperbolicity) ──────────────────────

def log_convexity_check(a, n_max):
    """
    d=2 hyperbolicity ↔ a_{n}² ≤ a_{n-1} × a_{n+1} (log-convexity when all a_n same sign).
    For alternating-sign GORZ_a_n (all positive): check a_n² ≤ a_{n-1} × a_{n+1}.
    """
    print("\nLog-convexity check: a_n^2 ≤ a_{n-1} × a_{n+1}?")
    violations = 0
    for n in range(1, min(n_max, len(a) - 1)):
        lhs = float(a[n])**2
        rhs = float(a[n-1]) * float(a[n+1])
        if lhs > rhs * (1 + 1e-10):  # allow tiny numerical error
            violations += 1
            if violations <= 3:
                print(f"  n={n}: a_{n}^2={lhs:.4e} > a_{n-1}×a_{n+1}={rhs:.4e}  VIOLATION")
    if violations == 0:
        print(f"  All log-convex for n=1..{min(n_max, len(a)-2)} ✓")
    else:
        print(f"  {violations} violations found ✗")
    return violations == 0

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("JENSEN POLYNOMIAL HYPERBOLICITY — NEWTON RECURRENCE")
    print("=" * 70)

    # Parameters
    N_ZEROS = 20000   # zeros to use for P_{2k}
    K_MAX   = 230     # max k for power sums (need k up to n_max)
    N_COEFFS = 225    # number of GORZ a_n coefficients to compute
    D_MAX   = 30      # max degree d to test
    N_MAX   = 200     # max n to test per d

    # Step 1: Load zeros
    print(f"\nStep 1: Load {N_ZEROS} zeros...")
    gammas = load_zeros(N_ZEROS)
    print(f"Loaded {len(gammas)} zeros. γ_1={gammas[0]:.6f}, γ_last={gammas[-1]:.6f}")

    # Step 2: Power sums
    P = compute_power_sums(gammas, K_MAX, add_tail=True)

    # Step 3: Newton recurrence for GORZ a_n
    a = compute_gorz_coefficients(P, N_COEFFS)

    # Step 4: Verify against known values
    print(f"\nStep 4: Verification")
    print(f"  a[0] = {float(a[0]):.10f}  (expected Ξ(0) ≈ 0.4971207782)")
    print(f"  a[1] = {float(a[1]):.10e}  (expected ≈ 1.14860e-02)")
    print(f"  a[2] = {float(a[2]):.10e}  (expected ≈ 1.23452e-04)")
    print(f"  Signs: {['+' if float(a[n])>0 else '-' for n in range(10)]}")

    # Step 5: Log-convexity
    log_convexity_check(a, min(N_COEFFS, 50))

    # Step 6: Discriminant for specific d values
    analyze_discriminant_growth(a, d=9, n_vals=[0, 5, 10, 20, 50, 100, 150, 200])

    # Step 7: Full hyperbolicity table
    table = run_hyperbolicity_table(a, D_MAX, N_MAX)

    # Summary
    print(f"\n{'='*70}")
    print(f"SUMMARY: Jensen Polynomial Hyperbolicity d=1..{D_MAX}, n=0..{N_MAX}")
    print(f"{'='*70}")
    max_d_all_hyp = max((d for d, r in table.items() if r["all_hyperbolic"]), default=0)
    min_d_fail = min((d for d, r in table.items() if not r["all_hyperbolic"]), default=None)
    print(f"  Max d with ALL n hyperbolic: {max_d_all_hyp}")
    print(f"  First d with any failure:     {min_d_fail}")
    if max_d_all_hyp >= 9:
        print(f"  ✓ GORZ extended from d≤8 to d≤{max_d_all_hyp} computationally!")
    else:
        print(f"  ✗ Could not extend GORZ beyond d=8")

    # Save
    out_path = os.path.join(RESULTS_DIR, "jensen_newton.json")
    summary = {
        "n_zeros": len(gammas),
        "k_max": K_MAX,
        "n_coeffs": N_COEFFS,
        "d_max": D_MAX,
        "n_max": N_MAX,
        "a0": float(a[0]),
        "a1": float(a[1]),
        "max_d_all_hyperbolic": max_d_all_hyp,
        "first_d_with_failure": min_d_fail,
        "table": {str(d): r for d, r in table.items()}
    }
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
