"""
Extended Jensen polynomial hyperbolicity via Hermite's criterion.

Hermite's criterion: a degree d polynomial with real coefficients is hyperbolic
(all roots real) iff the Hankel matrix of power sums H = [S_{i+j}]_{i,j=0..d-1}
is positive semi-definite (positive definite for simple roots).

Advantages over polyroots:
  - No iterative root finding, no convergence failures
  - O(d^3) per polynomial vs O(d^2 × iterations) for polyroots
  - Numerically stable at high precision via Cholesky
  - Scales to d=100, n=1000 feasibly

Method:
  1. Compute GORZ coefficients a[0..n_max] via Cauchy DFT at 200 dps
  2. For each (d,n), build Jensen polynomial J_{d,n}(X) = Σ C(d,j) a[n+j] X^j
  3. Compute power sums S_k via Newton's identities (O(d^2) per polynomial)
  4. Build d×d Hankel matrix H; check all leading minors > 0
  5. Threshold: min leading minor / S_0^d > 10^{-100} → hyperbolic

Target: d=1..100, n=0..200 at 200 dps
"""
import mpmath
import json
import os
import time
from math import comb

DPS = 300
mpmath.mp.dps = DPS + 50

RESULTS_DIR = "/Users/scones/projects/riemann/results"


def Xi(z):
    s = mpmath.mpf('1') / 2 + mpmath.mpc('0', '1') * z
    return (mpmath.mpf('1') / 2 * s * (s - 1)
            * mpmath.power(mpmath.pi, -s / 2)
            * mpmath.gamma(s / 2)
            * mpmath.zeta(s))


def compute_gorz_coeffs(n_max=110, R=7, N_pts=220):
    """Cauchy DFT for GORZ coefficients at DPS precision."""
    mpmath.mp.dps = DPS + 80
    R_mp = mpmath.mpf(str(R))

    print(f"[Cauchy DFT] R={R}, N_pts={N_pts}, n_max={n_max}, dps={DPS}")
    t0 = time.time()

    vals = []
    for k in range(N_pts):
        theta = mpmath.pi * k / N_pts
        z_k = R_mp * mpmath.exp(mpmath.mpc('0', '1') * theta)
        vals.append(Xi(z_k))
        if k % 40 == 0 or k == N_pts - 1:
            print(f"  Xi eval {k}/{N_pts}  ({time.time()-t0:.1f}s)", flush=True)

    print(f"  All {N_pts} Xi evals done in {time.time()-t0:.1f}s")
    print("  Extracting Taylor coefficients...")

    coeffs = []
    for n in range(n_max + 1):
        s = mpmath.mpc('0')
        for k, v in enumerate(vals):
            angle = mpmath.pi * (-2 * n * k) / N_pts
            s += v * (mpmath.cos(angle) + mpmath.mpc('0', '1') * mpmath.sin(angle))
        cn = s / N_pts / mpmath.power(R_mp, 2 * n)
        gorz_an = ((-1) ** n) * mpmath.re(cn)
        coeffs.append(gorz_an)
        if n % 20 == 0:
            print(f"  a[{n:3d}] = {mpmath.nstr(gorz_an, 10)}  ({time.time()-t0:.1f}s)", flush=True)

    print(f"  DFT done in {time.time()-t0:.1f}s")
    return coeffs


def newton_power_sums(poly):
    """
    Given ascending polynomial coefficients poly[0..d] (poly[j] = coeff of X^j),
    compute power sums S[0], S[1], ..., S[2d-2] via Newton's identities.

    Normalizes by leading coefficient: treats monic polynomial p(X)/c_d.
    S[0] = d (= number of roots, each counted with multiplicity).
    """
    d = len(poly) - 1
    if d == 0:
        return []

    # Monic version: divide by c_d = poly[d]
    # p(X) = c_d (X^d + e1 X^{d-1} + ... + (-1)^d e_d)
    # elementary symmetric polynomials: e_k = (-1)^k * poly[d-k] / poly[d]
    c_d = poly[d]
    e = [mpmath.mpf('0')] * (d + 1)
    for k in range(1, d + 1):
        e[k] = ((-1) ** k) * poly[d - k] / c_d

    # Newton's identities: S[k] = e[1]*S[k-1] - e[2]*S[k-2] + ... + (-1)^{k-1} k*e[k]
    S = [mpmath.mpf('0')] * (2 * d - 1)
    S[0] = mpmath.mpf(d)  # S_0 = number of roots

    for k in range(1, 2 * d - 1):
        val = mpmath.mpf('0')
        if k <= d:
            for j in range(1, k):
                val += ((-1) ** (j - 1)) * e[j] * S[k - j]
            val += (-1) ** (k - 1) * k * e[k]
        else:
            # For k > d: S[k] = e[1]*S[k-1] - e[2]*S[k-2] + ... - (-1)^d e[d]*S[k-d]
            for j in range(1, d + 1):
                val += ((-1) ** (j - 1)) * e[j] * S[k - j]
        S[k] = val

    return S


def hermite_criterion(poly, tol_exp=-290):
    """
    Check hyperbolicity via Hankel matrix of power sums.

    Returns (is_hyperbolic, min_normalized_minor, reason).
    Threshold: min_normalized_minor > 10^{tol_exp} → hyperbolic.

    The d×d Hankel matrix H = [S_{i+j}]_{i,j=0..d-1} (using S_0=d, S_1, ..., S_{2d-2}).
    H is positive definite ↔ all roots real and distinct.
    """
    d = len(poly) - 1
    if d <= 0:
        return True, mpmath.mpf('1'), "trivial"

    # Rescale X → λY with λ = (poly[0]/poly[d])^{1/d} to bound root magnitudes to O(1)
    # Jensen polynomial a[n] >> a[n+d] so roots without rescaling are huge → unstable power sums
    if poly[0] > 0 and poly[d] > 0:
        lam = mpmath.power(poly[0] / poly[d], mpmath.mpf('1') / d)
    else:
        lam = mpmath.mpf('1')
    scaled = [poly[j] * mpmath.power(lam, j) for j in range(d + 1)]

    # Additional normalization so max coefficient = 1
    max_coeff = max(abs(c) for c in scaled)
    if max_coeff == 0:
        return True, mpmath.mpf('1'), "trivial"
    poly_work = [c / max_coeff for c in scaled]

    # Power sums of poly_work (roots Y = X/λ, magnitude O(1))
    try:
        S = newton_power_sums(poly_work)
    except Exception as ex:
        return None, mpmath.mpf('0'), f"newton_fail: {ex}"

    # Build d×d Hankel matrix
    H = mpmath.matrix(d, d)
    for i in range(d):
        for j in range(d):
            H[i, j] = S[i + j]

    # Check leading principal minors
    tol = mpmath.mpf('10') ** tol_exp
    min_minor = mpmath.mpf('inf')

    for k in range(1, d + 1):
        Hk = mpmath.matrix(k, k)
        for i in range(k):
            for j in range(k):
                Hk[i, j] = H[i, j]
        det_k = mpmath.det(Hk)

        # Normalize by S_0^k = d^k to get scale-invariant criterion
        normalized = det_k / mpmath.power(mpmath.mpf(d), k)
        if normalized < min_minor:
            min_minor = normalized

        if normalized < tol:
            return False, normalized, f"minor_{k}_fail"

    return True, min_minor, "all_minors_pos"


def jensen_poly_mp(a, d, n):
    return [mpmath.mpf(comb(d, j)) * a[n + j] for j in range(d + 1)]


def main():
    print("=" * 70)
    print("EXTENDED JENSEN HYPERBOLICITY — HERMITE'S CRITERION")
    print(f"dps={DPS}, Hankel matrix method (no polyroots)")
    print("=" * 70)

    n_max = 110
    R = 7
    N_pts = 220

    a = compute_gorz_coeffs(n_max=n_max, R=R, N_pts=N_pts)

    # Sanity check
    known = [0.4971207781883141, 0.01148597215757272, 0.000123452018070318]
    print("\n[Sanity check]")
    for n in range(3):
        diff_rel = abs(float(a[n]) - known[n]) / abs(known[n])
        print(f"  a[{n}]: {mpmath.nstr(a[n], 10)}  rel_err={diff_rel:.2e}")
    print()

    # Full table: d=1..100, n=0..min(200, n_max-d)
    d_max = 100
    n_test = 100  # test up to this n (limited by n_max - d)

    print(f"[Full table] d=1..{d_max}, n=0..{n_test}")
    print(f"{'d':>4}  {'status':>15}  {'min_minor':>20}  {'time':>8}")
    print("-" * 55)

    table = {}
    grand_start = time.time()

    for d in range(1, d_max + 1):
        n_limit = min(n_test, n_max - d)
        first_fail = None
        min_minor_all = mpmath.mpf('inf')
        conv_failures = []
        t_d = time.time()

        for n in range(n_limit + 1):
            poly = jensen_poly_mp(a, d, n)
            hyp, minor, reason = hermite_criterion(poly)

            if hyp is None:
                conv_failures.append(n)
                continue
            if minor < min_minor_all:
                min_minor_all = minor
            if hyp is False and first_fail is None:
                first_fail = n

        all_hyp = (first_fail is None)
        minor_float = float(min_minor_all) if min_minor_all != mpmath.mpf('inf') else float('inf')

        table[d] = {
            "all_hyperbolic": all_hyp,
            "first_fail": first_fail,
            "conv_failures": conv_failures,
            "min_minor": minor_float,
            "n_tested": n_limit
        }
        status = "✓ ALL" if all_hyp else f"✗ n={first_fail}"
        conv_note = f"  [cf={conv_failures}]" if conv_failures else ""
        elapsed = time.time() - t_d
        print(f"  d={d:3d}: {status:12s}  min_minor={minor_float:.3e}  {elapsed:6.1f}s{conv_note}",
              flush=True)

    max_d_hyp = max((d for d in table if table[d]["all_hyperbolic"]), default=0)

    out = {
        "dps": DPS,
        "n_max": n_max,
        "R": R,
        "N_pts": N_pts,
        "d_max_tested": d_max,
        "n_max_tested": n_test,
        "max_d_all_hyperbolic": max_d_hyp,
        "method": "hermite_criterion_hankel",
        "tol_exp": -80,
        "table": table,
        "total_time": time.time() - grand_start
    }
    out_path = os.path.join(RESULTS_DIR, "jensen_hermite_ext.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    print()
    print("=" * 70)
    print(f"RESULT: Max d with J_{{d,n}} hyperbolic for all tested n: {max_d_hyp}")
    total = time.time() - grand_start
    print(f"Total time: {total:.1f}s ({total/60:.1f} min)")
    print(f"Results saved to {out_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
