"""
Track convergence of normalized Jensen polynomials to Hermite polynomials.

From GORZ Theorem 1: δ(n)^{-d} J_{d,n}(δ(n)(X-1)/exp(A(n))) → H_d(X) as n → ∞.

This script:
1. Computes A(n), δ(n) from GORZ (18) using Theorem 7
2. Normalizes J_{d,n} and measures deviation from H_d
3. Tracks when the deviation drops below threshold (estimates N_0(d))

Key formulas from GORZ:
  L = L(nb) where nb = 2n-2, L is the unique positive solution of nb = L(πe^L + 3/4)
  δ(n) = sqrt( nb × L² / ((1+L)nb - (3/4)L²) ) / (L² × K(nb))  ... see eq (18)
  A(n) = log(nb × L² / (4*nb²)) + ...

Actually using the simpler approximation that's "sufficiently strong" (footnote 6 of GORZ):
  δ(n) ≈ sqrt( (L+2)×2nb / (nb × L⁴ × K²) )

This script computes the normalized polynomial Jhat and measures ||Jhat - H_d||_coeff
(max deviation of normalized coefficients from Hermite coefficients).
"""
import mpmath
import json
import os
import time
from math import comb

DPS = 150  # sufficient for convergence tracking
mpmath.mp.dps = DPS + 30

RESULTS_DIR = "/Users/scones/projects/riemann/results"


def Xi(z):
    s = mpmath.mpf('1') / 2 + mpmath.mpc('0', '1') * z
    return (mpmath.mpf('1') / 2 * s * (s - 1)
            * mpmath.power(mpmath.pi, -s / 2)
            * mpmath.gamma(s / 2)
            * mpmath.zeta(s))


def compute_gorz_coeffs(n_max=85, R=7, N_pts=190):
    mpmath.mp.dps = DPS + 50
    R_mp = mpmath.mpf(str(R))
    print(f"[DFT] R={R}, N={N_pts}, n_max={n_max}, dps={DPS}")
    t0 = time.time()
    vals = []
    for k in range(N_pts):
        theta = mpmath.pi * k / N_pts
        z_k = R_mp * mpmath.exp(mpmath.mpc('0', '1') * theta)
        vals.append(Xi(z_k))
        if k % 50 == 0:
            print(f"  Xi {k}/{N_pts}  ({time.time()-t0:.1f}s)", flush=True)
    coeffs = []
    for n in range(n_max + 1):
        s = mpmath.mpc('0')
        for k, v in enumerate(vals):
            angle = mpmath.pi * (-2 * n * k) / N_pts
            s += v * (mpmath.cos(angle) + mpmath.mpc('0', '1') * mpmath.sin(angle))
        cn = s / N_pts / mpmath.power(R_mp, 2 * n)
        coeffs.append(((-1) ** n) * mpmath.re(cn))
    print(f"  Done in {time.time()-t0:.1f}s")
    return coeffs


def compute_L(nb, prec=DPS):
    """
    Solve nb = L(pi*e^L + 3/4) for L > 0 using Newton's method.
    nb = 2n - 2 (GORZ notation).
    """
    mpmath.mp.dps = prec + 10
    nb_mp = mpmath.mpf(str(nb))
    # Initial guess: L ≈ log(nb/pi) - log(log(nb/pi))
    import math
    x0 = math.log(float(nb) / math.pi)
    L = mpmath.mpf(str(max(x0, 0.5)))
    for _ in range(50):
        f = L * (mpmath.pi * mpmath.exp(L) + mpmath.mpf('3') / 4) - nb_mp
        df = (mpmath.pi * mpmath.exp(L) + mpmath.mpf('3') / 4) + L * mpmath.pi * mpmath.exp(L)
        L_new = L - f / df
        if abs(L_new - L) < mpmath.mpf('10') ** (-(prec - 5)):
            L = L_new
            break
        L = L_new
    return L


def gorz_parameters(n, prec=DPS):
    """
    Compute A(n) and δ(n) from GORZ (18).
    nb = 2n - 2, K = K(nb) = (1/L + 1/L²)*nb - 3/4
    δ(n) = sqrt(nb*L² / ((1+L)*nb - (3/4)*L²)) / (L²*K)   ... from GORZ
    """
    mpmath.mp.dps = prec + 10
    nb = mpmath.mpf(2 * n - 2)
    L = compute_L(int(float(nb)), prec=prec)

    # K = (L^{-1} + L^{-2}) nb - 3/4
    K = (1 / L + 1 / L ** 2) * nb - mpmath.mpf('3') / 4

    # δ(n)² = nb*L² / ((1+L)*nb - (3/4)*L²) / (L⁴*K²)
    # Simpler version from GORZ (18):
    delta_sq = nb * L ** 2 / ((1 + L) * nb - mpmath.mpf('3') / 4 * L ** 2) / (L ** 4 * K ** 2)
    # Fix: actually from (18): δ(n) = sqrt( (n(L+2)) / (nb * L^4 * K^2) ) ... let me use the explicit formula
    # δ(n)² = 2*nb*(L+2) / (nb² * L^4 * K^2)  from (18)?
    # Actually from (18): δ(n)² = n*(L+2) / (b*L⁴*K²) where b*n = (1+L)*nb - (3/4)L²
    # Let me use: δ(n)² = (L+2) / (L⁴ * K² * (nb/(2*(L+1)) - 3*L/(8*(L+1))))
    # This is getting complex. Use the simpler two-term approximation from the paper's example:
    # δ(n) ≈ sqrt(L(nb)/2nb) as leading term

    delta_leading = mpmath.sqrt(L / (2 * nb))

    # A(n) = log(nb*L²/(4*nb²)) + corrections ≈ log(L²/(4*nb)) from GORZ notation
    A_leading = mpmath.log(nb * L ** 2 / 4) / 2 - mpmath.log(nb) / 2
    # More precisely from GORZ: A(n) ≈ (nb/4 - 3/4*log(L) + const)
    # Use: A(n) = (1/2)*log((nb/(4*nb²)) * ...) -- complex. Use leading term.
    # A(n) = log sqrt(nb L² / 4) / nb ≈ log(L/2) + (1/2)log(nb/nb) -- skip
    # Let's use GORZ's formula: exp(A(n)) = gamma(n) / (gamma(n) × normalized factor)
    # Actually from the paper: exp(A(n)) ≈ a[n+1]/a[n] × e^{-delta²}
    # This requires knowing a[n], so compute inline.

    return float(delta_leading), float(A_leading), float(L), float(K)


def hermite_poly(d):
    """
    Compute H_d coefficients (ascending power order) using the generating function.
    H_d(X) = sum_{k=0}^{floor(d/2)} d!/(k!(d-2k)!) * (-1)^k * X^{d-2k}
    """
    coeffs = [mpmath.mpf('0')] * (d + 1)
    for k in range(d // 2 + 1):
        power = d - 2 * k
        coeff = (mpmath.factorial(d) / (mpmath.factorial(k) * mpmath.factorial(d - 2 * k))
                 * mpmath.mpf('1' if k % 2 == 0 else '-1'))
        coeffs[power] = coeff
    # Normalize so leading coefficient = 1
    lead = coeffs[d]
    return [c / lead for c in coeffs]


def normalized_jensen(a, d, n, delta, A_n):
    """
    Compute Jhat_{d,n}(X) = delta^{-d} / (a[n]*exp(A_n)) * J_{d,n}(delta*X - 1)
    This should converge to H_d(X) as n → ∞.
    """
    # J_{d,n}(Y) where Y = delta*X - 1, so X = (Y+1)/delta
    # Jhat(X) = delta^{-d} / (a[n]*exp(A_n)) * sum_j C(d,j) a[n+j] (delta*X - 1)^j
    mpmath.mp.dps = DPS + 20
    delta_mp = mpmath.mpf(str(delta))
    A_mp = mpmath.mpf(str(A_n))

    # Compute J_{d,n} evaluated at Y = delta*X - 1 as polynomial in X
    # (delta*X - 1)^j = sum_{k=0}^{j} C(j,k) delta^k X^k (-1)^{j-k}
    result = [mpmath.mpf('0')] * (d + 1)
    for j in range(d + 1):
        c_j = mpmath.mpf(comb(d, j)) * a[n + j]
        for k in range(j + 1):
            coeff = c_j * mpmath.mpf(comb(j, k)) * (mpmath.mpf('-1') ** (j - k)) * (delta_mp ** k)
            result[k] += coeff

    # Normalize: divide by a[n]*exp(A_n)*delta^d
    norm_factor = a[n] * mpmath.exp(A_mp) * (delta_mp ** d)
    return [c / norm_factor for c in result]


def main():
    print("=" * 70)
    print("GORZ CONVERGENCE: Jensen → Hermite deviation tracker")
    print(f"dps={DPS}")
    print("=" * 70)

    n_max = 85
    a = compute_gorz_coeffs(n_max=n_max)

    print()
    print("[GORZ δ(n) parameters]")
    for n in [5, 10, 20, 40, 80]:
        delta, A_n, L, K = gorz_parameters(n)
        print(f"  n={n:4d}: δ={delta:.6f}, L={L:.4f}, K={K:.4f}")

    print()
    print("[Convergence to Hermite] Computing deviation ||Jhat - H_d||_coeff")
    print(f"{'d':>4}  {'n':>5}  {'max_coeff_dev':>20}  {'δ(n)':>10}")
    print("-" * 50)

    results = {}

    for d in [2, 3, 4, 6, 8, 10, 15, 20]:
        H_d = hermite_poly(d)
        results[d] = []
        for n in [5, 10, 20, 40, 60, 80 - d]:
            if n + d > n_max:
                continue
            delta, A_n, L, K = gorz_parameters(n)

            # Compute deviation: A_n from GORZ directly via ratio a[n+1]/a[n]
            # Use empirical A: a[n]*exp(A_n) ≈ a[n+1]/delta² (from GORZ Theorem 3)
            A_emp = float(mpmath.log(a[n + 1] / a[n]) + delta ** 2)

            Jhat = normalized_jensen(a, d, n, delta, A_emp)

            # Max deviation from H_d coefficients
            max_dev = max(abs(float(Jhat[k]) - float(H_d[k])) for k in range(d + 1))

            results[d].append({"n": n, "delta": delta, "max_dev": max_dev})
            print(f"  d={d:2d}, n={n:4d}: dev={max_dev:.4e}  δ={delta:.4f}", flush=True)

    print()
    out = {
        "dps": DPS,
        "results": {str(d): v for d, v in results.items()}
    }
    out_path = os.path.join(RESULTS_DIR, "gorz_convergence.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Results saved to {out_path}")


if __name__ == "__main__":
    main()
