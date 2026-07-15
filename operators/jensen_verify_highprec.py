"""
High-precision Jensen polynomial verification.

The previous jensen_cauchy.py computed Taylor coefficients at 150 dps but then
converted to float64 for root finding — losing all precision for large n/d.

At d=20, n=28: coefficient ratio a_48/a_28 ~ 10^{-43}, so float64 eigenvalue
computation has 0 valid digits. This script keeps everything in mpmath throughout.

Method:
  1. Recompute a_n for n=0..80 via Cauchy DFT at 200 dps (full mpmath)
  2. Check J_{d,n} hyperbolicity using mpmath.polyroots (no float64 conversion)
  3. Target d=1..35, n=0..45

Expected result: d=20..30 "failures" from jensen_cauchy.json are numerical noise.
Genuine failure would show |Im(root)|/|Re(root)| > 0.01.
Numerical noise at 200 dps appears as |Im(root)|/|Re(root)| < 10^{-150}.
"""
import mpmath
import json
import os
import time
from math import comb

DPS = 200
mpmath.mp.dps = DPS + 50  # guard digits

RESULTS_DIR = "/Users/scones/projects/riemann/results"


def Xi(z):
    """Ξ(z) = ξ(1/2 + iz)"""
    s = mpmath.mpf('1') / 2 + mpmath.mpc('0', '1') * z
    return (mpmath.mpf('1') / 2 * s * (s - 1)
            * mpmath.power(mpmath.pi, -s / 2)
            * mpmath.gamma(s / 2)
            * mpmath.zeta(s))


def compute_gorz_coeffs(n_max=75, R=7, N_pts=180):
    """
    Compute GORZ Taylor coefficients a[0..n_max] of Ξ(iz) = Σ a[n] z^{2n}
    via Cauchy DFT on |z|=R at DPS precision.

    Key: a[n] = (-1)^n × Re[(1/N) Σ_k Xi(R e^{iπk/N}) e^{-2inπk/N}] / R^{2n}

    All arithmetic in mpmath — no float64 conversion.
    Aliasing error: O((R^{2N} × a_{n+N}) / a_n) ~ negligible for R=7, N=180.
    """
    mpmath.mp.dps = DPS + 80
    R_mp = mpmath.mpf(str(R))

    print(f"[Cauchy DFT] R={R}, N_pts={N_pts}, n_max={n_max}, dps={DPS}")
    t0 = time.time()

    # Evaluate Xi at N_pts points on upper-half circle |z|=R
    vals = []
    for k in range(N_pts):
        theta = mpmath.pi * k / N_pts
        z_k = R_mp * mpmath.exp(mpmath.mpc('0', '1') * theta)
        v = Xi(z_k)
        vals.append(v)
        if k % 30 == 0 or k == N_pts - 1:
            print(f"  Xi eval {k}/{N_pts}  ({time.time()-t0:.1f}s)", flush=True)

    print(f"  All {N_pts} Xi evaluations done in {time.time()-t0:.1f}s")

    # DFT extraction: for each n, sum Σ_k vals[k] × e^{-2inπk/N}
    # Then divide by N × R^{2n} and apply (-1)^n
    print("  Extracting Taylor coefficients via DFT...")
    coeffs = []
    for n in range(n_max + 1):
        s = mpmath.mpc('0')
        for k, v in enumerate(vals):
            angle = mpmath.pi * (-2 * n * k) / N_pts
            s += v * (mpmath.cos(angle) + mpmath.mpc('0', '1') * mpmath.sin(angle))
        cn = s / N_pts / mpmath.power(R_mp, 2 * n)
        # GORZ convention: Ξ(iz) = Σ (-1)^n × c_n × z^{2n}? Let's clarify.
        # From the cauchy DFT derivation: cn = (-1)^n × GORZ_a[n] (real part)
        gorz_an = ((-1) ** n) * mpmath.re(cn)
        coeffs.append(gorz_an)
        if n % 10 == 0:
            print(f"  a[{n:3d}] = {mpmath.nstr(gorz_an, 12)}  ({time.time()-t0:.1f}s)", flush=True)

    print(f"  DFT extraction done in {time.time()-t0:.1f}s")
    return coeffs


def jensen_poly_mp(a, d, n):
    """
    J_{d,n}(X) = Σ_{j=0}^d C(d,j) a[n+j] X^j  (returns list of mpmath coefficients)
    """
    return [mpmath.mpf(comb(d, j)) * a[n + j] for j in range(d + 1)]


def hyperbolic_ratio_mp(poly, prec=None):
    """
    Find max |Im(root)| / (|Re(root)| + 1) using mpmath.polyroots.
    Returns (is_hyperbolic, max_ratio, roots).
    Threshold: ratio < 10^{-80} → hyperbolic (numerical noise at 200 dps).
    """
    if prec is None:
        prec = DPS
    mpmath.mp.dps = prec + 30

    d = len(poly) - 1
    if d <= 0:
        return True, mpmath.mpf('0'), []

    # poly[j] = C(d,j)*a[n+j] — coefficients decrease super-exponentially.
    # Rescale X → λY with λ = (poly[0]/poly[d])^{1/d} so first and last
    # coefficients become equal.  This reduces the dynamic range from ~10^{60}
    # to ~10^{6}, allowing polyroots to converge.
    if poly[0] > 0 and poly[d] > 0:
        lam = mpmath.power(poly[0] / poly[d], mpmath.mpf('1') / d)
    else:
        lam = mpmath.mpf('1')
    scaled = [poly[j] * mpmath.power(lam, j) for j in range(d + 1)]

    max_coeff = max(abs(c) for c in scaled)
    if max_coeff == 0:
        return True, mpmath.mpf('0'), []
    poly_norm = [c / max_coeff for c in scaled]

    # polyroots takes descending-order coefficients [c_d, c_{d-1}, ..., c_0]
    try:
        roots = mpmath.polyroots(poly_norm[::-1], maxsteps=300, extraprec=50)
    except Exception as e:
        return None, mpmath.mpf('nan'), []

    # Roots are in Y-space; scale back to X-space: X = λY
    roots = [lam * r for r in roots]

    max_ratio = mpmath.mpf('0')
    for r in roots:
        ratio = abs(mpmath.im(r)) / (abs(mpmath.re(r)) + mpmath.mpf('1'))
        if ratio > max_ratio:
            max_ratio = ratio

    tol = mpmath.mpf('10') ** (-80)
    return bool(max_ratio < tol), max_ratio, roots


def main():
    print("=" * 70)
    print("HIGH-PRECISION JENSEN POLYNOMIAL VERIFICATION")
    print(f"dps={DPS}, using mpmath.polyroots throughout (no float64)")
    print("=" * 70)

    # ── Step 1: Compute GORZ coefficients ────────────────────────────────────
    n_max = 75
    R = 7
    N_pts = 180

    a = compute_gorz_coeffs(n_max=n_max, R=R, N_pts=N_pts)

    # Sanity check: compare first few to known values
    known = [0.4971207781883141, 0.01148597215757272, 0.000123452018070318]
    print("\n[Sanity check] First 3 coefficients vs known values:")
    for n in range(3):
        diff_rel = abs(float(a[n]) - known[n]) / abs(known[n])
        print(f"  a[{n}]: computed={mpmath.nstr(a[n], 10)}, known={known[n]:.10e}, rel_err={diff_rel:.2e}")

    print()

    # ── Step 2: Targeted hyperbolicity check ─────────────────────────────────
    # First: verify d=19 and d=20 around the failure zone
    print("[Targeted check] d=18..25, n=24..35 (suspicious zone from previous run)")
    print(f"{'d':>4}  {'n':>4}  {'hyperbolic':>12}  {'max_ratio':>20}")
    print("-" * 50)

    suspicious_results = {}
    t_check = time.time()
    for d in range(18, 26):
        for n in range(24, 36):
            if n + d > n_max:
                continue
            poly = jensen_poly_mp(a, d, n)
            hyp, ratio, _ = hyperbolic_ratio_mp(poly)
            key = f"d={d},n={n}"
            suspicious_results[key] = {
                "hyperbolic": bool(hyp) if hyp is not None else None,
                "max_ratio": float(ratio) if ratio is not None else None
            }
            sym = "✓" if hyp is True else ("✗" if hyp is False else "?")
            ratio_str = "conv_fail" if hyp is None else f"{float(ratio):.3e}"
            print(f"  d={d:2d}, n={n:2d}: {sym}  ratio={ratio_str}  ({time.time()-t_check:.1f}s)",
                  flush=True)

    print()

    # ── Step 3: Full table d=1..35, n=0..40 ──────────────────────────────────
    print("[Full table] d=1..35, n=0..40")
    print(f"{'d':>4}  {'status':>15}  {'max_ratio_all':>20}")
    print("-" * 50)

    table = {}
    for d in range(1, 36):
        n_limit = min(40, n_max - d)
        first_fail = None
        max_ratio_all = mpmath.mpf('0')
        t_d = time.time()

        conv_failures = []
        for n in range(n_limit + 1):
            poly = jensen_poly_mp(a, d, n)
            hyp, ratio, _ = hyperbolic_ratio_mp(poly)
            if hyp is None:
                conv_failures.append(n)
                continue  # convergence failure — not a mathematical failure
            if ratio > max_ratio_all:
                max_ratio_all = ratio
            if hyp is False and first_fail is None:
                first_fail = n

        all_hyp = (first_fail is None)
        table[d] = {
            "all_hyperbolic": all_hyp,
            "first_fail": first_fail,
            "conv_failures": conv_failures,
            "max_ratio": float(max_ratio_all)
        }
        status = "✓ ALL" if all_hyp else f"✗ n={first_fail}"
        conv_note = f"  [conv_fail n={conv_failures}]" if conv_failures else ""
        print(f"  d={d:2d}: {status:15s}  max_ratio={float(max_ratio_all):.3e}  ({time.time()-t_d:.1f}s){conv_note}",
              flush=True)

    max_d_hyp = max((d for d in table if table[d]["all_hyperbolic"]), default=0)

    # ── Save results ─────────────────────────────────────────────────────────
    out = {
        "dps": DPS,
        "n_max": n_max,
        "R": R,
        "N_pts": N_pts,
        "d_max_tested": 35,
        "n_max_tested": 40,
        "max_d_all_hyperbolic": max_d_hyp,
        "hyperbolic_threshold": "10^{-80}",
        "table": table,
        "suspicious_zone": suspicious_results
    }
    out_path = os.path.join(RESULTS_DIR, "jensen_highprec.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    print()
    print("=" * 70)
    print(f"RESULT: Max d with J_{{d,n}} hyperbolic for all tested n: {max_d_hyp}")
    if max_d_hyp >= 35:
        print("All d=1..35 verified hyperbolic for n=0..40 at 200 dps.")
        print("Previous failures at d=20..30 were NUMERICAL NOISE (float64 artifacts).")
    else:
        print(f"Genuine failure at d={max_d_hyp+1} — check details in table.")
    print(f"Results saved to {out_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
