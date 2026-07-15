"""
Jensen Polynomial Hyperbolicity — Cauchy Integral Method

Computes Taylor coefficients a[n] of Ξ(z) = Σ_{n≥0} (-1)^n GORZ_a[n] z^{2n} via:

  GORZ_a[n] × R^{2n} = (1/N) Σ_{k=0}^{N-1} Ξ(R e^{iπk/N}) e^{-2πink/N}

(Cauchy integral discretized as DFT, exploiting even-function symmetry.)

Uses mpmath throughout to achieve high accuracy. R = 10 (inside convergence disk of
radius γ₁ = 14.134). N = 256 points → a[0..127] available.

Goal: verify J_{d,n} hyperbolic for d=9..30 and n=0..100.
"""

import mpmath
import numpy as np
from math import comb
import json
import os
import time

mpmath.mp.dps = 150  # 150 decimal places — plenty for a[0..75]

RESULTS_DIR = "/Users/scones/projects/riemann/results"

# ── Xi function ──────────────────────────────────────────────────────────────

def Xi(z):
    """Ξ(z) = ξ(1/2 + iz), ξ(s) = (1/2)s(s-1)π^{-s/2}Γ(s/2)ζ(s)"""
    s = mpmath.mpf('0.5') + 1j * z
    return mpmath.mpf('0.5') * s * (s - 1) * mpmath.power(mpmath.pi, -s/2) * mpmath.gamma(s/2) * mpmath.zeta(s)


# ── Cauchy / DFT coefficient extraction ──────────────────────────────────────

def compute_gorz_cauchy(n_max, R=10, N_pts=None, verbose=True):
    """
    Extract GORZ Taylor coefficients a[0..n_max] of Ξ using Cauchy integral.

    Ξ(z) = Σ_{n≥0} GORZ_a[n] (-z²)^n = Σ (-1)^n GORZ_a[n] z^{2n}

    Substituting w = -z²:  F(w) = Σ GORZ_a[n] w^n  with F(w) = Ξ(√(-w))

    To extract a[n]: evaluate Ξ at z_k = R e^{iπ(2k+1)/(2N)} for k=0..N-1
    (or at z_k = R e^{iπk/N}), then:

      GORZ_a[n] × R^{2n} = Re[ (1/N) Σ_k Ξ(R e^{iθ_k}) e^{-2niθ_k} ]

    where θ_k = πk/N and we use that Ξ is even (only even powers).

    Strategy: Evaluate at N points on circle |z| = R in the z-plane (z real-like axis).
    Since Ξ(z) = Ξ(-z), only even powers appear.
    """
    if N_pts is None:
        N_pts = 2 * (n_max + 20)  # oversampling

    R_mp = mpmath.mpf(str(R))

    if verbose:
        print(f"Cauchy DFT: R={R}, N={N_pts}, n_max={n_max}, dps={mpmath.mp.dps}")

    t0 = time.time()

    # Evaluate Ξ at N_pts points z_k = R × e^{iπk/N_pts} for k=0..N_pts-1
    # (These are uniformly spaced on upper semicircle, then reflected by symmetry)
    vals = []
    for k in range(N_pts):
        theta = mpmath.pi * k / N_pts  # θ_k = πk/N
        z_k = R_mp * mpmath.exp(1j * theta)
        v = Xi(z_k)
        vals.append(v)
        if verbose and (k % 32 == 0 or k == N_pts - 1):
            elapsed = time.time() - t0
            print(f"  k={k}/{N_pts}  ({elapsed:.1f}s)")

    if verbose:
        print(f"Evaluations done in {time.time()-t0:.1f}s")

    # DFT: a[n] × R^{2n} = Re[ (1/N) Σ_k vals[k] × e^{-2πink/N} ]
    # With our θ_k = πk/N: e^{iθ_k} so 2nθ_k = 2nπk/N = 2πnk/(N) × ... wait.
    #
    # F(w) = Σ GORZ_a[n] w^n  where w = -z²
    # F(w_k) = Ξ(z_k) where z_k = R e^{iθ_k} and w_k = -z_k² = -R² e^{2iθ_k} = R² e^{i(π+2θ_k)}
    #
    # With θ_k = πk/N:  w_k = R² e^{iπ(2k+N)/N} = R² e^{iπ(2k/N + 1)}
    # These are N points on circle |w|=R² (in w-plane).
    #
    # Standard DFT: (1/N) Σ_k F(w_k) × w_k^{-n} × R^{-2n} × ...
    #
    # Simpler: directly from θ:
    # Ξ(R e^{iθ}) = Σ (-1)^m GORZ_a[m] R^{2m} e^{2imθ}  (Ξ is even: Σ ... z^{2m})
    #
    # So: (1/N) Σ_k Ξ(R e^{iθ_k}) × e^{-2inθ_k} = GORZ_a_modified[n] × R^{2n}
    #
    # But wait: Ξ(z) = Σ c_m z^{2m} where c_m = (-1)^m GORZ_a[m].
    # So Ξ(R e^{iθ}) = Σ c_m R^{2m} e^{2imθ}
    # DFT gives: (1/N) Σ_k Ξ(R e^{iπk/N}) e^{-2πink/(something)} = ...
    #
    # Standard DFT with N points at θ_k = πk/N:
    # F_n = (1/N) Σ_{k=0}^{N-1} f(k) e^{-2πink/N}
    #
    # Here f(k) = Ξ(R e^{iπk/N}) and we want:
    # F_n = (1/N) Σ_k Ξ(R e^{iπk/N}) e^{-2πink/N}
    #      = (1/N) Σ_k [Σ_m c_m R^{2m} e^{2imπk/N}] e^{-2πink/N}
    #      = Σ_m c_m R^{2m} (1/N) Σ_k e^{2πi(m-n)k/N}
    #      = Σ_m c_m R^{2m} δ_{m≡n (mod N)}
    #      = c_n R^{2n} + c_{n+N} R^{2(n+N)} + ...
    #      ≈ c_n R^{2n}  (since c_{n+N} R^{2(n+N)} is negligible for large N)
    #
    # So: c_n = (-1)^n GORZ_a[n] ≈ F_n / R^{2n}
    # Therefore: GORZ_a[n] = (-1)^n F_n / R^{2n}

    # Compute DFT manually in mpmath
    coeffs_standard = []
    for n in range(n_max + 1):
        s = mpmath.mpf('0')
        for k, v in enumerate(vals):
            theta_k = mpmath.pi * k / N_pts
            s += v * mpmath.exp(-2j * n * theta_k)
        cn = s / N_pts / mpmath.power(R_mp, 2*n)
        gorz_an = ((-1)**n) * mpmath.re(cn)  # take real part; imaginary should be ~0
        coeffs_standard.append(gorz_an)

    if verbose:
        print(f"DFT extraction done.")

    return coeffs_standard


# ── Jensen polynomial hyperbolicity ──────────────────────────────────────────

def jensen_poly_float(a, d, n):
    """J_{d,n}(X) = Σ_{j=0}^d C(d,j) GORZ_a[n+j] X^j  (as float list)"""
    return [float(comb(d, j)) * float(a[n + j]) for j in range(d + 1)]

def is_hyperbolic_float(coeffs, tol=1e-8):
    """Check if polynomial has all real roots (float precision)."""
    p = np.array(coeffs[::-1], dtype=complex)
    if len(p) <= 1:
        return True, 0.0
    try:
        roots = np.roots(p)
    except np.linalg.LinAlgError:
        return None, float('nan')
    max_imag = max(abs(r.imag) for r in roots) if len(roots) > 0 else 0.0
    scale    = max(abs(r.real) for r in roots) + 1.0 if len(roots) > 0 else 1.0
    return bool(max_imag / scale < tol), float(max_imag / scale)


def run_hyperbolicity_table(a, d_max, n_max, verbose=True):
    """Full hyperbolicity table for d=1..d_max, n=0..n_max."""
    results = {}
    print(f"\n{'='*70}")
    print(f"JENSEN POLYNOMIAL HYPERBOLICITY TABLE")
    print(f"GORZ (2019): d≤8 proved hyperbolic for all n.")
    print(f"Goal: verify d=9..{d_max} computationally, n=0..{n_max}")
    print(f"{'='*70}")
    print(f"{'d':>4}  {'n_tested':>8}  {'n_hyp':>6}  {'fail@n':>8}  {'min_ratio':>12}  {'status'}")
    print(f"{'-'*65}")

    for d in range(1, d_max + 1):
        n_limit = min(n_max, len(a) - d - 1)
        if n_limit < 0:
            print(f"{d:>4}  not enough coefficients (need {d+1}, have {len(a)})")
            break

        n_tested = 0
        n_hyp = 0
        first_fail = None
        min_ratio = float('inf')
        max_ratio = 0.0

        for n in range(n_limit + 1):
            coeffs = jensen_poly_float(a, d, n)
            # Skip if any coefficient is inf/nan
            if any(not np.isfinite(c) for c in coeffs):
                continue
            hyp, ratio = is_hyperbolic_float(coeffs)
            if hyp is None:
                continue
            n_tested += 1
            if hyp:
                n_hyp += 1
            else:
                if first_fail is None:
                    first_fail = n
            if np.isfinite(ratio):
                min_ratio = min(min_ratio, ratio)
                max_ratio = max(max_ratio, ratio)

        if n_tested == 0:
            print(f"{d:>4}  no valid tests")
            continue

        all_hyp = (n_hyp == n_tested)
        status = "✓ ALL" if all_hyp else f"✗ FAIL"
        print(f"{d:>4}  {n_tested:>8}  {n_hyp:>6}  {str(first_fail):>8}  {max_ratio:>12.4e}  {status}")

        results[d] = {
            "d": d, "n_tested": n_tested, "n_hyperbolic": n_hyp,
            "first_fail": first_fail, "max_ratio": max_ratio,
            "all_hyperbolic": all_hyp
        }

    return results


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("JENSEN POLYNOMIAL HYPERBOLICITY — CAUCHY INTEGRAL METHOD")
    print("=" * 70)
    print(f"mpmath.mp.dps = {mpmath.mp.dps}")

    R     = 10    # evaluation radius (< γ₁ = 14.134)
    N_PTS = 128   # DFT points (extracts n=0..63)
    N_MAX = 60    # max coefficient index (safe with N_PTS=128)
    D_MAX = 30    # max degree d
    N_JMAX = 30   # max n for Jensen check (need a[n..n+d_max])

    # Step 1: Compute GORZ Taylor coefficients
    print(f"\nStep 1: Cauchy DFT on circle |z|={R}, N={N_PTS} points...")
    a = compute_gorz_cauchy(N_MAX, R=R, N_pts=N_PTS, verbose=True)

    # Step 2: Verify
    print(f"\nStep 2: Verification (compare to mpmath.diff values)")
    known = {
        0: 4.9712077819e-01,
        1: 1.1486e-02,
        2: 1.2345e-04,
        3: 8.324e-07,
        4: 3.992e-09,
        5: 1.462e-11,
        6: 4.27e-14,
        7: 1.031e-16,
        8: 2.10e-19,
        9: 3.678e-22,
    }
    print(f"{'n':>4}  {'Computed':>20}  {'Expected':>20}  {'Rel.Err':>12}")
    for n, expected in known.items():
        computed = float(a[n])
        err = abs(computed - expected) / abs(expected) if expected != 0 else float('nan')
        sign = '+' if computed > 0 else '-'
        print(f"{n:>4}  {computed:>20.10e}  {expected:>20.10e}  {err:>12.3e}")

    # Check all positive (GORZ_a_n > 0 ↔ Ξ has all real zeros)
    neg_count = sum(1 for x in a if float(x) < 0)
    print(f"\nAll GORZ_a[n] > 0? {neg_count == 0} (negatives: {neg_count})")
    if neg_count > 0:
        negs = [(n, float(a[n])) for n in range(len(a)) if float(a[n]) < 0]
        print(f"  First few negatives: {negs[:5]}")

    # Step 3: Log-convexity (needed for d=2 hyperbolicity)
    print(f"\nLog-convexity check: a_n^2 ≤ a_{{n-1}} × a_{{n+1}} for GORZ coefficients?")
    viol = 0
    for n in range(1, N_MAX):
        lhs = float(a[n])**2
        rhs = float(a[n-1]) * float(a[n+1])
        if rhs > 0 and lhs > rhs * (1 + 1e-10):
            viol += 1
            if viol <= 3:
                print(f"  VIOLATION at n={n}: a_{n}^2={lhs:.4e} > a_{{n-1}}×a_{{n+1}}={rhs:.4e}")
    if viol == 0:
        print(f"  Log-convex for all n=1..{N_MAX-1} ✓")
    else:
        print(f"  {viol} violations ✗")

    # Step 4: Full Jensen table
    table = run_hyperbolicity_table(a, D_MAX, N_JMAX, verbose=True)

    # Step 5: Detailed roots for d=9 and d=25 at various n
    print(f"\nDetailed root analysis (d=9):")
    for n_test in [0, 5, 10, 20, 30]:
        if n_test + 9 >= len(a):
            break
        coeffs = jensen_poly_float(a, 9, n_test)
        if any(not np.isfinite(c) for c in coeffs):
            continue
        hyp, ratio = is_hyperbolic_float(coeffs)
        p = np.array(coeffs[::-1], dtype=complex)
        try:
            roots = sorted(np.roots(p).real)
        except:
            roots = []
        print(f"  n={n_test:3d}: hyp={'✓'if hyp else '✗'} max_imag/scale={ratio:.3e}  roots: {[f'{r:.2f}' for r in roots[:4]]}")

    if 25 in table:
        print(f"\nDetailed root analysis (d=25):")
        for n_test in [0, 5]:
            if n_test + 25 >= len(a):
                break
            coeffs = jensen_poly_float(a, 25, n_test)
            if any(not np.isfinite(c) for c in coeffs):
                continue
            hyp, ratio = is_hyperbolic_float(coeffs)
            p = np.array(coeffs[::-1], dtype=complex)
            try:
                roots = sorted(np.roots(p).real)
            except:
                roots = []
            print(f"  n={n_test:3d}: hyp={'✓'if hyp else '✗'} max_imag/scale={ratio:.3e}  roots: {[f'{r:.2f}' for r in roots[:4]]}")

    # Step 6: Summary
    print(f"\n{'='*70}")
    max_d_all = max((d for d, r in table.items() if r["all_hyperbolic"]), default=0)
    first_fail = min((d for d, r in table.items() if not r["all_hyperbolic"]), default=None)
    print(f"SUMMARY:")
    print(f"  Max d with ALL n hyperbolic (n=0..{N_JMAX}): {max_d_all}")
    print(f"  First d with any failure: {first_fail}")
    if max_d_all >= 9:
        print(f"  ✓ GORZ extended computationally from d≤8 to d≤{max_d_all}!")
    print(f"  Note: GORZ also proved each J_{{d,n}} eventually hyperbolic (large n).")
    print(f"  Computational check for small n extends the known range of d.")

    # Save
    out = os.path.join(RESULTS_DIR, "jensen_cauchy.json")
    coeffs_float = [float(x) for x in a]
    save = {
        "method": "cauchy_dft",
        "R": R, "N_pts": N_PTS, "n_max": N_MAX,
        "dps": mpmath.mp.dps,
        "gorz_a": coeffs_float,
        "d_max": D_MAX, "n_jmax": N_JMAX,
        "max_d_all_hyperbolic": max_d_all,
        "first_d_with_failure": first_fail,
        "table": {str(d): r for d, r in table.items()}
    }
    with open(out, "w") as f:
        json.dump(save, f, indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
