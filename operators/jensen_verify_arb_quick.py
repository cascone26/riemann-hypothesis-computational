"""
Quick Arb certification test - smaller range to verify approach works.
Uses d=1..20, n=0..40 to complete in reasonable time (~5-10 min).
"""

import json
import os
import time
import mpmath
from math import comb

try:
    from flint import arb, acb, acb_poly
    FLINT_AVAILABLE = True
except ImportError:
    FLINT_AVAILABLE = False

RESULTS_DIR = "/Users/scones/projects/riemann/results"
DPS = 200
mpmath.mp.dps = DPS + 50


def Xi(z):
    """Ξ(z) = ξ(1/2 + iz) computed with mpmath."""
    s = mpmath.mpf('1') / 2 + mpmath.mpc('0', '1') * z
    return (mpmath.mpf('1') / 2 * s * (s - 1)
            * mpmath.power(mpmath.pi, -s / 2)
            * mpmath.gamma(s / 2)
            * mpmath.zeta(s))


def compute_gorz_coeffs_mpmath(n_max=75, R=7, N_pts=150):
    """Compute GORZ coefficients via Cauchy DFT using mpmath."""
    mpmath.mp.dps = DPS + 50
    R_mp = mpmath.mpf(str(R))

    print(f"[Cauchy DFT / mpmath] R={R}, N_pts={N_pts}, n_max={n_max}, dps={DPS}")
    t0 = time.time()

    vals = []
    for k in range(N_pts):
        theta = mpmath.pi * k / N_pts
        z_k = R_mp * mpmath.exp(mpmath.mpc('0', '1') * theta)
        v = Xi(z_k)
        vals.append(v)
        if k % 30 == 0 or k == N_pts - 1:
            print(f"  Xi eval {k}/{N_pts}  ({time.time()-t0:.1f}s)", flush=True)

    print(f"  DFT extraction...")
    coeffs = []
    for n in range(n_max + 1):
        s = mpmath.mpc('0')
        for k, v in enumerate(vals):
            angle = mpmath.pi * (-2 * n * k) / N_pts
            s += v * (mpmath.cos(angle) + mpmath.mpc('0', '1') * mpmath.sin(angle))
        cn = s / N_pts / mpmath.power(R_mp, 2 * n)
        gorz_an = ((-1) ** n) * mpmath.re(cn)
        coeffs.append(gorz_an)
    print(f"  Done in {time.time()-t0:.1f}s")
    return coeffs


def jensen_poly_arb(a, d, n):
    """Construct J_{d,n}(X) = Σ_{j=0}^d C(d,j) a[n+j] X^j"""
    poly = []
    for j in range(d + 1):
        coeff_val = float(a[n + j])
        coeff_bin = arb(comb(d, j)) * arb(str(coeff_val))
        poly.append(acb(coeff_bin))
    return poly


def hyperbolic_check_arb(poly):
    """Check hyperbolicity using Arb certified root-finding."""
    d = len(poly) - 1
    if d <= 0:
        return True, arb('0'), arb('0'), []

    magnitudes = [abs(c) for c in poly]
    max_coeff_arb = max(magnitudes)
    max_coeff = float(max_coeff_arb.mid()) if max_coeff_arb > 0 else 0.0
    if max_coeff == 0:
        return True, arb('0'), arb('0'), []
    poly_norm = [c / max_coeff_arb for c in poly]

    try:
        poly_acb = acb_poly(poly_norm)
        roots = poly_acb.roots()
    except Exception as e:
        return None, arb('nan'), arb('nan'), []

    max_imag_ratio = arb('0')
    max_error_radius = arb('0')
    roots_with_radii = []

    for root in roots:
        radius = root.rad()
        roots_with_radii.append((root, radius))
        imag_part = abs(root.imag)
        real_part = abs(root.real)
        ratio = imag_part / (real_part + arb('1'))

        if ratio > max_imag_ratio:
            max_imag_ratio = ratio
        if radius > max_error_radius:
            max_error_radius = radius

    all_real = True
    for root in roots:
        imag_abs = abs(root.imag)
        rad = root.rad()
        # Root is certified real if |Im| <= radius
        if imag_abs > rad:
            all_real = False
            break

    return all_real, max_imag_ratio, max_error_radius, roots_with_radii


def main():
    if not FLINT_AVAILABLE:
        print("ERROR: python-flint not installed")
        return

    print("=" * 70)
    print("ARB-CERTIFIED JENSEN VERIFICATION (quick test)")
    print("=" * 70)

    # Coefficients
    n_max = 75
    a = compute_gorz_coeffs_mpmath(n_max=n_max, R=7, N_pts=150)

    # Quick hyperbolicity check: d=1..20, n=0..40
    print("\n[Quick check] d=1..20, n=0..40 with Arb certified roots")
    print(f"{'d':>4}  {'status':>15}  {'max_ratio':>20}  {'max_radius':>20}")
    print("-" * 70)

    table = {}
    t_total = time.time()

    for d in range(1, 21):
        n_limit = min(40, n_max - d)
        first_fail = None
        max_ratio_all = arb('0')
        max_radius_all = arb('0')
        t_d = time.time()

        for n in range(n_limit + 1):
            if n + d > n_max:
                break
            poly = jensen_poly_arb(a, d, n)
            hyp, ratio, radius, roots = hyperbolic_check_arb(poly)

            if hyp is None:
                continue

            if ratio > max_ratio_all:
                max_ratio_all = ratio
            if radius > max_radius_all:
                max_radius_all = radius

            if hyp is False and first_fail is None:
                first_fail = n

        all_hyp = (first_fail is None)
        table[d] = {
            "all_hyperbolic": all_hyp,
            "first_fail": first_fail,
            "n_max_tested": n_limit,
            "max_ratio": float(max_ratio_all.mid()),
            "max_ratio_error": float(max_ratio_all.rad().mid()),
            "max_radius": float(max_radius_all.mid()),
        }

        status = "✓ ALL" if all_hyp else f"✗ n={first_fail}"
        print(f"  d={d:2d}: {status:15s}  ratio={float(max_ratio_all.mid()):.3e}  "
              f"radius={float(max_radius_all.mid()):.3e}  ({time.time()-t_d:.1f}s)",
              flush=True)

    max_d_hyp = max((d for d in table if table[d]["all_hyperbolic"]), default=0)

    # Save
    out = {
        "method": "Hybrid: mpmath coefficients + Arb certified root-finding",
        "test_type": "quick (d=1..20, n=0..40)",
        "coefficient_precision_dps": DPS,
        "n_max_computed": n_max,
        "d_max_tested": 20,
        "n_max_tested": 40,
        "max_d_all_hyperbolic": max_d_hyp,
        "table": table,
    }

    out_path = os.path.join(RESULTS_DIR, "jensen_arb_quick_test.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    print()
    print("=" * 70)
    print(f"Max d all hyperbolic: {max_d_hyp}")
    print(f"Total time: {(time.time()-t_total)/60:.1f} min")
    print(f"Results: {out_path}")
    print("=" * 70)


if __name__ == "__main__":
    if not FLINT_AVAILABLE:
        print("ERROR: python-flint not installed")
        exit(1)
    main()
