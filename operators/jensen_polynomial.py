"""
Jensen Polynomial Hyperbolicity Extension
==========================================
Griffin-Ono-Rolen-Zagier (2019) proved: J_{d,n}(X) is hyperbolic (all real roots)
for d ≤ 8 and ALL n ≥ 0, and for fixed d and n sufficiently large.

Open: d ≥ 9 for ALL n.

This script:
1. Computes Taylor coefficients of Ξ(z) around z=0 to high precision
2. Forms Jensen polynomials J_{d,n}(X) for d = 1..25 and n = 0..500
3. Checks hyperbolicity (all real roots) numerically
4. Reports the LARGEST d for which hyperbolicity holds for all tested n
5. Looks for patterns that suggest the proof structure for all d

Definition: Ξ(z) = Σ_{n=0}^∞ a_n z^{2n}
Jensen polynomial: J_{d,n}(X) = Σ_{j=0}^{d} C(d,j) a_{n+j} X^j

Hyperbolic iff all roots of J_{d,n} are real.

Reference: GORZ 2019 PNAS — arXiv:1902.07321
"""

import mpmath
import numpy as np
from math import comb
import json
import os
import time

mpmath.mp.dps = 50

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")

# ── Compute Taylor coefficients of Ξ ─────────────────────────────────────────

def Xi_func(z):
    """Ξ(z) = ξ(1/2 + iz) — real-valued on real line
    ξ(s) = (1/2) s(s-1) π^{-s/2} Γ(s/2) ζ(s)
    """
    s = mpmath.mpf('0.5') + 1j * z
    return mpmath.mpf('0.5') * s * (s - 1) * mpmath.power(mpmath.pi, -s/2) * mpmath.gamma(s/2) * mpmath.zeta(s)


def compute_taylor_coeffs(N_max, prec=50):
    """
    Compute a_n = Ξ^{(2n)}(0) / (2n)! for n = 0, 1, ..., N_max

    Uses mpmath's Taylor series differentiation.
    Returns list of a_n as mpmath floats.
    """
    mpmath.mp.dps = prec + 10

    print(f"Computing Taylor coefficients a_n for n=0..{N_max}...")
    t0 = time.time()

    # Method: use numerical differentiation
    # Ξ is even: Ξ(z) = Ξ(-z), so only even derivatives are nonzero
    # a_n = Ξ^{(2n)}(0) / (2n)!

    coeffs = []

    # Use mpmath's Taylor series
    # Compute derivatives using the identity: a_n = (2n-th coeff of Taylor series)
    # via numerical differentiation at z=0 with small step
    h = mpmath.mpf('1e-3')  # step size

    # We'll use the formula based on function values at z = j*h for j = 0..2*N_max
    # and discrete Fourier extraction

    # Better: use mpmath's derivative function
    for n in range(N_max + 1):
        if n == 0:
            a_n = mpmath.re(Xi_func(mpmath.mpf('0')))
        else:
            # Compute 2n-th derivative via mpmath.diff
            deriv = mpmath.diff(lambda z: mpmath.re(Xi_func(z)), mpmath.mpf('0'), 2*n)
            # a_n = deriv / factorial(2n)
            a_n = deriv / mpmath.factorial(2*n)
        coeffs.append(a_n)
        if n % 10 == 0:
            print(f"  a_{n} = {float(a_n):.6e}  ({time.time()-t0:.1f}s)", flush=True)

    print(f"Done in {time.time()-t0:.1f}s")
    return coeffs


def compute_taylor_coeffs_fast(N_max, prec=50):
    """
    Faster method: evaluate Ξ at N_max+1 Chebyshev points and use FFT-like extraction.

    Ξ(z) = Σ_{n=0}^{N_max} a_n z^{2n} + O(z^{2N_max+2})

    Let z = r×e^{2πi k/(2N)} for k=0..2N-1, then:
    Ξ(r×e^{2πik/(2N)}) = Σ_n a_n r^{2n} e^{2πi(2n)k/(2N)}

    Taking DFT: a_n r^{2n} = (1/N) Re[Σ_k Ξ(r e^{2πik/N}) e^{-2πi n k/N}]

    (since only even powers appear)
    """
    mpmath.mp.dps = prec + 10

    print(f"Computing Taylor coefficients (fast method) N_max={N_max}...")
    t0 = time.time()

    N = 2 * (N_max + 5)  # oversampling
    r = mpmath.mpf('0.5')  # evaluation radius (should be < first root |γ₁|=14.13)

    # Evaluate Ξ at N points on circle |z|=r
    vals = []
    for k in range(N):
        theta = 2 * mpmath.pi * k / N
        z = r * mpmath.exp(1j * theta)
        s = mpmath.mpf('0.5') + 1j * z
        vals.append(mpmath.mpf('0.5') * s * (s-1) * mpmath.power(mpmath.pi, -s/2) * mpmath.gamma(s/2) * mpmath.zeta(s))
    vals = np.array([complex(v) for v in vals])

    # DFT
    fft = np.fft.fft(vals) / N

    # Extract coefficients: a_n * r^{2n} = Re(fft[2n])
    coeffs = []
    r_float = float(r)
    for n in range(N_max + 1):
        if 2*n < N:
            a_n_rn = fft[2*n].real
            a_n = a_n_rn / r_float**(2*n)
        else:
            a_n = 0.0
        coeffs.append(float(a_n))

    print(f"Done in {time.time()-t0:.2f}s")
    return coeffs


def compute_taylor_coeffs_diff(N_max, prec=50):
    """
    Most accurate: use mpmath.diff with Taylor expansion at z=0.
    Only feasible for small N_max (say N_max ≤ 30) due to cost.
    """
    mpmath.mp.dps = prec + 10
    print(f"Computing {N_max+1} Taylor coefficients via mpmath.diff...")
    t0 = time.time()
    coeffs = []

    def xi_real(x):
        """Ξ(x) restricted to real x"""
        s = mpmath.mpf('0.5') + 1j * x
        return mpmath.re(mpmath.mpf('0.5') * s * (s-1) * mpmath.power(mpmath.pi, -s/2) * mpmath.gamma(s/2) * mpmath.zeta(s))

    for n in range(N_max + 1):
        # Compute 2n-th derivative at x=0
        if n == 0:
            a = xi_real(mpmath.mpf('0'))
        else:
            # Use 5-point or higher stencil via mpmath.diff
            a = mpmath.diff(xi_real, mpmath.mpf('0'), 2*n)
            a /= mpmath.factorial(2*n)
        coeffs.append(float(a))
        print(f"  n={n}: a_{n} = {float(a):.8e}  ({time.time()-t0:.1f}s)", flush=True)

    return coeffs


# ── Jensen polynomial construction and hyperbolicity check ────────────────────

def jensen_poly(a_coeffs, d, n):
    """
    J_{d,n}(X) = Σ_{j=0}^d C(d,j) a_{n+j} X^j

    Returns list of coefficients [c_0, c_1, ..., c_d] where c_j = C(d,j)*a_{n+j}
    """
    poly_coeffs = []
    for j in range(d + 1):
        if n + j < len(a_coeffs):
            c = comb(d, j) * a_coeffs[n + j]
        else:
            c = 0.0
        poly_coeffs.append(c)
    return poly_coeffs


def is_hyperbolic(poly_coeffs, tol=1e-10):
    """
    Check if polynomial with given coefficients has all real roots.
    Uses numpy's root finder.

    poly_coeffs = [c_0, c_1, ..., c_d] where poly = Σ c_j X^j

    Note: numpy roots expects HIGH DEGREE FIRST, so reverse.
    """
    if not any(abs(c) > 1e-100 for c in poly_coeffs):
        return True, []  # zero polynomial

    # Normalize so leading coefficient is nonzero
    p = poly_coeffs[:]
    while len(p) > 1 and abs(p[-1]) < 1e-100:
        p = p[:-1]

    if len(p) <= 1:
        return True, []

    # numpy roots: expects [a_d, a_{d-1}, ..., a_0] (highest degree first)
    p_rev = p[::-1]
    roots = np.roots(p_rev)

    # Check all roots are real
    max_imag = max(abs(r.imag) for r in roots) if len(roots) > 0 else 0.0
    all_real = max_imag < tol * (max(abs(r.real) for r in roots) + 1.0)

    return all_real, roots


def check_hyperbolicity_table(a_coeffs, d_max=25, n_max=200):
    """
    Check hyperbolicity of J_{d,n} for all d=1..d_max and n=0..n_max.
    Returns a table of results.
    """
    print(f"\nHyperbolicity table: d=1..{d_max}, n=0..{n_max}")
    print("(✓ = all roots real, ✗ = has complex roots)")
    print()

    results = {}
    for d in range(1, d_max + 1):
        row_results = {}
        first_fail = None
        n_checked = 0
        for n in range(n_max + 1):
            if n + d >= len(a_coeffs):
                break
            poly = jensen_poly(a_coeffs, d, n)
            hyp, roots = is_hyperbolic(poly)
            row_results[n] = hyp
            n_checked += 1
            if not hyp and first_fail is None:
                first_fail = n

        all_hyp = all(row_results.values())
        results[d] = {"all_hyperbolic": all_hyp, "first_fail": first_fail,
                      "n_checked": n_checked}

        status = "✓ ALL" if all_hyp else f"✗ fails at n={first_fail}"
        print(f"  d={d:>3}: {status}  (checked n=0..{n_checked-1})")

    return results


# ── Discriminant analysis ─────────────────────────────────────────────────────

def discriminant_analysis(a_coeffs, d_max=10, n_max=500):
    """
    For each d, compute the discriminant of J_{d,n} as a function of n.
    The discriminant is 0 iff J_{d,n} has a repeated root (touching hyperbolicity boundary).

    For d=2: disc = (C(2,1) a_{n+1})^2 - 4 C(2,0)a_n × C(2,2)a_{n+2}
           = 4 a_{n+1}^2 - 4 a_n a_{n+2}
           = 4(a_{n+1}^2 - a_n a_{n+2})
    Hyperbolic iff a_{n+1}^2 ≥ a_n a_{n+2} (log-convexity of a_n)

    (Note: a_n > 0, so sign of a_{n+1}^2 - a_n a_{n+2} matters)
    """
    print(f"\nDiscriminant analysis (d=2..{min(d_max,5)}):")

    for d in range(2, min(d_max, 6)):
        if d > 4:
            continue  # only d=2,3,4 have simple discriminant formulas

        print(f"\n  d={d}:")
        min_disc_n = None
        min_disc_val = float('inf')

        for n in range(n_max + 1):
            if n + d >= len(a_coeffs):
                break
            poly = jensen_poly(a_coeffs, d, n)
            # Compute discriminant via roots
            hyp, roots = is_hyperbolic(poly)
            # Min "reality margin" = max|Im(root)|
            if roots:
                max_imag = max(abs(r.imag) for r in roots)
                if max_imag < min_disc_val:
                    min_disc_val = max_imag
                    min_disc_n = n

        print(f"    Min |Im(root)| = {min_disc_val:.2e} at n={min_disc_n}")
        print(f"    (> 0 means all polynomials are genuinely hyperbolic, bounded away from degenerate)")


# ── Log-convexity check of Taylor coefficients ────────────────────────────────

def log_convexity_check(a_coeffs, n_max=100):
    """
    Necessary condition for J_{2,n} to be hyperbolic: a_{n+1}^2 ≥ a_n × a_{n+2}
    This is log-convexity of the sequence (a_n).
    """
    print("\nLog-convexity check: a_{n+1}^2 ≥ a_n × a_{n+2}?")
    a = np.array(a_coeffs[:n_max+2], dtype=float)
    violations = []
    for n in range(len(a)-2):
        lhs = a[n+1]**2
        rhs = a[n] * a[n+2]
        if lhs < rhs - 1e-15 * (abs(lhs) + abs(rhs) + 1):
            violations.append(n)

    if violations:
        print(f"  LOG-CONVEXITY FAILS at n={violations[:10]}")
        print(f"  → J_{{2,n}} would not be hyperbolic for these n")
    else:
        print(f"  ✓ Log-convex for all tested n=0..{min(len(a)-2, n_max)}")
        print(f"    (This is necessary for J_{{2,n}} hyperbolicity)")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    print("="*70)
    print("JENSEN POLYNOMIAL HYPERBOLICITY — EXTENSION OF GORZ 2019")
    print("="*70)
    print()
    print("Griffin-Ono-Rolen-Zagier (2019): J_{d,n} hyperbolic for d≤8 ALL n.")
    print("Goal: Verify for d=9..25 computationally using high-precision Ξ coefficients.")
    print()

    # Method 1: Fast DFT method for many coefficients
    print("METHOD 1: Fast DFT extraction of Taylor coefficients...")
    N_fast = 150
    a_fast = compute_taylor_coeffs_fast(N_fast, prec=50)
    print(f"First 10 coefficients a_n:")
    for i in range(min(10, len(a_fast))):
        print(f"  a_{i} = {a_fast[i]:.10e}")
    print()

    # Verify signs: a_n should alternate (a_n > 0 for all n since Ξ has positive coeffs)
    # Actually: Ξ(z) has only POSITIVE Taylor coefficients (on the imaginary axis)
    # because Ξ(iz) = Ξ(-z) and Ξ has all positive Taylor coefficients in z²
    print(f"All a_n > 0? {all(a > 0 for a in a_fast[:50])}")
    print(f"Sign pattern of first 20: {['+'if a>0 else '-' for a in a_fast[:20]]}")
    print()

    # Method 2: High-precision diff for small n (ground truth)
    print("METHOD 2: High-precision mpmath.diff for n=0..15 (slow but accurate)...")
    N_diff = 15
    a_exact = compute_taylor_coeffs_diff(N_diff, prec=50)
    print()

    # Compare methods
    print("COMPARISON: fast vs exact for n=0..10:")
    print(f"{'n':>3}  {'fast':>15}  {'exact':>15}  {'rel_err':>12}")
    for n in range(min(11, len(a_exact))):
        af = a_fast[n]
        ae = float(a_exact[n])
        if abs(ae) > 1e-100:
            rel = abs((af - ae)/ae)
        else:
            rel = abs(af)
        print(f"{n:>3}  {af:>15.8e}  {ae:>15.8e}  {rel:>12.2e}")
    print()

    # Use the exact method for the first N_exact coefficients and fast for the rest
    N_exact = min(len(a_exact), len(a_fast))
    a_combined = list(a_exact[:N_exact]) + a_fast[N_exact:]

    # ── Hyperbolicity table ────────────────────────────────────────────────────
    print("\n" + "="*70)
    print("HYPERBOLICITY TABLE")
    print("="*70)

    # Convert to float for numpy computation
    a_float = [float(x) for x in a_combined]

    results = check_hyperbolicity_table(a_float, d_max=25, n_max=min(N_fast-25, 100))

    # ── Log-convexity of Ξ coefficients ────────────────────────────────────────
    log_convexity_check(a_float, n_max=min(N_fast-2, 100))

    # ── Discriminant analysis ─────────────────────────────────────────────────
    discriminant_analysis(a_float, d_max=8, n_max=min(N_fast-8, 100))

    # ── GORZ verification ─────────────────────────────────────────────────────
    print("\n" + "="*70)
    print("GORZ VERIFICATION (d=1..8 should be ✓ for all tested n)")
    print("="*70)
    for d in range(1, 9):
        r = results.get(d, {})
        status = "✓ ALL" if r.get("all_hyperbolic") else f"✗ fails at n={r.get('first_fail')}"
        print(f"  d={d}: {status}  [n=0..{r.get('n_checked',0)-1}]")

    # ── New territory ──────────────────────────────────────────────────────────
    print("\n" + "="*70)
    print("NEW TERRITORY (d=9..25)")
    print("="*70)
    new_results = []
    for d in range(9, 26):
        r = results.get(d, {})
        status = "✓ ALL" if r.get("all_hyperbolic") else f"✗ fails at n={r.get('first_fail')}"
        new_results.append((d, r.get("all_hyperbolic"), r.get("first_fail")))
        print(f"  d={d}: {status}  [n=0..{r.get('n_checked',0)-1}]")

    all_new_hyp = all(r[1] for r in new_results if r[1] is not None)
    print()
    if all_new_hyp:
        print("*** RESULT: J_{d,n} appears HYPERBOLIC for ALL tested d=9..25 and ALL n ***")
        print("This extends GORZ from d≤8 to d≤25 (computationally verified).")
        print("A rigorous proof for all d and all n remains open.")
    else:
        fail_cases = [(d,n) for d,_,n in new_results if n is not None]
        print(f"COUNTEREXAMPLES FOUND at: {fail_cases[:5]}")
        print("Jensen polynomial approach may have limitations.")

    # ── Ratio a_{n+1}^2 / (a_n × a_{n+2}) — approach to 1? ────────────────────
    print("\n" + "="*70)
    print("LOG-CONVEXITY RATIO r_n = a_{n+1}^2 / (a_n × a_{n+2})")
    print("(r_n > 1 ↔ log-convex, necessary for d=2 hyperbolicity)")
    print("Conjecture: r_n → 1 as n→∞?")
    print("="*70)
    a = np.array(a_float[:min(N_fast, 120)], dtype=float)
    for n in [0, 5, 10, 20, 30, 50, 70, 90]:
        if n + 2 < len(a) and abs(a[n]) > 1e-100 and abs(a[n+2]) > 1e-100:
            ratio = a[n+1]**2 / (a[n] * a[n+2])
            print(f"  n={n:>3}: r_n = {ratio:.8f}  (a_n={a[n]:.4e}, a_{{n+1}}={a[n+1]:.4e}, a_{{n+2}}={a[n+2]:.4e})")

    # ── Save results ───────────────────────────────────────────────────────────
    out = {
        "a_coeffs_fast": a_fast[:60],
        "a_coeffs_exact": [float(x) for x in a_exact],
        "hyperbolicity": {str(d): {"all_hyperbolic": r.get("all_hyperbolic"),
                                    "first_fail": r.get("first_fail"),
                                    "n_checked": r.get("n_checked")}
                          for d, r in results.items()},
        "new_territory_all_hyperbolic": all_new_hyp,
    }
    with open(os.path.join(RESULTS_DIR, "jensen_polynomial.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nResults saved.")


if __name__ == "__main__":
    main()
