#!/usr/bin/env python3
"""
gaussian_family_positivity.py — SuperPRISM 2026-07-20, item: Gaussian test-function family.

NEW, non-circular test of the single target inequality (research/uniform_positivity_target_2026_05_19.md):

    (*)  A(g) + B(g) >= P(g)   for every admissible even g   <=>  RH

Never tried in this project: pick h(r) := ghat(1/2+ir) = exp(-r^2/w^2) (Gaussian in the FREQUENCY
variable). Its position-space Fourier dual is ALSO Gaussian:

    g_w(x) = (w/(2*sqrt(pi))) * exp(-w^2 x^2 / 4)

This is the classical explicit-formula test family (Poitou / Odlyzko-te Riele style smoothing), never
run before against OUR non-circular tau = A+B-P construction.

Why it matters: w sweeps a ONE-PARAMETER family interpolating between the two regimes LEDGER §34c
identified as the crux:
  - w -> 0 :  h concentrates on r in [-O(w),O(w)] -- probes ONLY the sign-indefinite Omega core
              (|r|<6.29). g_w spreads out in x -- P(g) sums MANY primes (slowly decaying).
  - w -> oo:  h spreads over all r (including the coercive Omega>0 tail). g_w -> delta at x=0 --
              P(g) collapses to the first few primes.
Finite-N bandlimited-basis numerics (lambda^2=14) never isolated this axis directly; this closed-form
family sweeps it exactly, all the way into the core-dominated regime.

TERMS:
  A(g_w) = (1/2pi) int Omega(r) exp(-r^2/w^2) dr        [robust adaptive quadrature; Omega zero-free]
  B(g_w) = h(i/2) + h(-i/2) = 2*exp(1/(4 w^2))            [poles of xi at s=0,1 -- classical]
  P(g_w) = 2 * sum_{k>=2} Lambda(k) k^{-1/2} g_w(log k)   [von Mangoldt, prime powers]

Cross-check: A also computed via the Mittag-Leffler/Voigt closed form
  Omega(r) = C0 + sum_n [1/(n+1) - P_n(r)], P_n(r)=a_n/(a_n^2+r^2/4), a_n=n+1/4
  (1/2pi) int P_n(r) exp(-r^2/w^2) dr = erfcx(2 a_n / w)
against the quadrature result at several widths (verify-before-ship).
"""

import numpy as np
from scipy.integrate import quad
from scipy.special import digamma, erfcx
import json
import time

GAMMA_E = 0.5772156649015329
LOG_PI = float(np.log(np.pi))
C0 = -GAMMA_E - LOG_PI
SQRT_PI = float(np.sqrt(np.pi))


def omega(r):
    return np.real(digamma(0.25 + 0.5j * r)) - LOG_PI


def A_gauss_quad(w):
    """A(g_w) via adaptive quadrature. Range sized so exp(-r^2/w^2) has decayed to ~1e-16."""
    R = w * 6.1 + 5.0  # +5 floor so tiny w still resolves the O(1) core structure
    integrand = lambda r: omega(r) * np.exp(-r * r / (w * w))
    val, err = quad(integrand, -R, R, limit=400, epsabs=1e-12, epsrel=1e-10)
    return val / (2 * np.pi), err / (2 * np.pi)


def A_gauss_voigt_crosscheck(w, n_terms=20000):
    """Independent closed-form cross-check via the Mittag-Leffler/Voigt (erfcx) series."""
    n = np.arange(n_terms)
    a_n = n + 0.25
    term = (w / (2 * SQRT_PI * (n + 1))) - erfcx(2 * a_n / w)
    total = np.sum(term)
    # leading asymptotic tail correction: term_n ~ -0.75*w/(2*sqrt(pi)) / n^2 for n >> w
    tail_correction = -0.75 * w / (2 * SQRT_PI) * (1.0 / n_terms)
    return C0 * w / (2 * SQRT_PI) + total + tail_correction


def B_gauss(w):
    """B(g_w) = 2*exp(1/(4w^2)), the rank-2 pole term at s=0,1."""
    return 2.0 * np.exp(1.0 / (4.0 * w * w))


def primes_upto(n):
    """Fast numpy sieve of Eratosthenes."""
    if n < 2:
        return np.array([], dtype=np.int64)
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    return np.nonzero(sieve)[0].astype(np.int64)


# hard feasibility cap on prime enumeration (50M) — below w~0.7 the natural cutoff exceeds this
KMAX_HARD = 50_000_000
_PRIME_CACHE = {"kmax": 0, "primes": None, "logp": None}


def get_primes(kmax):
    kmax = min(kmax, KMAX_HARD)
    if _PRIME_CACHE["kmax"] < kmax:
        pr = primes_upto(kmax)
        _PRIME_CACHE.update(kmax=kmax, primes=pr, logp=np.log(pr.astype(np.float64)))
    m = _PRIME_CACHE["primes"] <= kmax
    return _PRIME_CACHE["primes"][m], _PRIME_CACHE["logp"][m]


def kmax_natural(w, floor_log=1e-16):
    """Natural prime cutoff where g_w(log k) drops below floor_log (uncapped)."""
    bound = (2.0 / w) * np.sqrt(-np.log(floor_log))
    return np.exp(bound)


def P_gauss(w):
    """P(g_w) = 2 sum_k Lambda(k) k^{-1/2} g_w(log k) over prime powers k=p^m.
    Returns (P, truncated_flag, natural_kmax)."""
    nat = kmax_natural(w)
    kmax = int(min(nat, KMAX_HARD))
    kmax = max(kmax, 50)
    truncated = nat > KMAX_HARD
    primes, logp = get_primes(kmax)
    coef = w / (2 * SQRT_PI)

    total = 0.0
    # m=1 (primes) vectorized
    x = logp
    total += np.sum(logp * (primes.astype(np.float64) ** -0.5) * coef * np.exp(-w * w * x * x / 4.0))
    # higher prime powers p^m <= kmax, m>=2 (sparse; iterate)
    pf = primes.astype(np.float64)
    m = 2
    while True:
        pm = pf ** m
        mask = pm <= kmax
        if not np.any(mask):
            break
        xm = m * logp[mask]
        total += np.sum(logp[mask] * (pm[mask] ** -0.5) * coef * np.exp(-w * w * xm * xm / 4.0))
        m += 1
    return 2.0 * total, truncated, nat


def main():
    print("=" * 78)
    print("GAUSSIAN-FAMILY WEIL POSITIVITY SCAN (non-circular; A by quadrature, cross-checked)")
    print("=" * 78)

    # --- cross-check quadrature vs Voigt closed form at a few widths first ---
    print("\n-- cross-check: quadrature vs Mittag-Leffler/Voigt closed form --")
    for w_test in [0.1, 1.0, 10.0, 100.0]:
        Aq, err = A_gauss_quad(w_test)
        Av = A_gauss_voigt_crosscheck(w_test)
        print(f"  w={w_test:7.2f}  A_quad={Aq:14.6f} (err~{err:.1e})  A_voigt={Av:14.6f}  "
              f"diff={Aq - Av:+.2e}")
    print()

    # feasible range: brute-force prime side needs ~e^(12/w) primes; hard cap 50M => w >~ 0.66.
    # widths below ~0.7 are prime-tail-truncated (flagged, not silently dropped).
    widths = sorted(set(
        list(np.geomspace(0.30, 2.0, 22)) +
        list(np.geomspace(2.0, 50.0, 16)) +
        list(np.geomspace(50.0, 2000.0, 10))
    ))

    results = []
    t0 = time.time()
    n_trunc = 0

    for w in widths:
        Aw, Aerr = A_gauss_quad(w)
        Bw = B_gauss(w)
        Pw, truncated, nat = P_gauss(w)
        tau_w = Aw + Bw - Pw
        if truncated:
            n_trunc += 1

        rec = {
            "w": w, "A": Aw, "A_quad_err": Aerr, "B": Bw, "P": Pw, "tau": tau_w,
            "prime_truncated": bool(truncated), "natural_kmax": float(nat),
        }
        results.append(rec)
        flag = " [P TRUNCATED — under-resolved]" if truncated else ""
        print(f"  w={w:9.4f}  A={Aw:14.6f}  B={Bw:14.6f}  P={Pw:14.6f}  "
              f"tau={tau_w:+14.6f}{flag}", flush=True)

    dt = time.time() - t0
    # only trust widths where the prime side is fully resolved
    clean = [r for r in results if not r["prime_truncated"]]
    taus_clean = np.array([r["tau"] for r in clean])
    ws_clean = np.array([r["w"] for r in clean])
    min_idx = int(np.argmin(taus_clean))

    print("\n" + "=" * 78)
    print(f"Scanned {len(widths)} widths in {dt:.1f}s. {n_trunc} had truncated prime side (w too small).")
    print(f"Over the {len(clean)} FULLY-RESOLVED widths (prime side exact):")
    print(f"  MIN tau(w) = {taus_clean[min_idx]:+.6f} at w = {ws_clean[min_idx]:.4f}")
    print(f"  Sign: {'ALL POSITIVE (no violation)' if taus_clean.min() > 0 else '*** NEGATIVE tau FOUND ***'}")
    print(f"Smallest fully-resolved width reached: w = {ws_clean.min():.4f} "
          f"(core |r|<6.29 probed for r-content up to ~{ws_clean.min()*6:.1f})")
    print("=" * 78)

    out = {
        "results": results,
        "min_tau": float(taus_clean.min()),
        "min_tau_w": float(ws_clean[min_idx]),
        "all_positive": bool(taus_clean.min() > 0),
        "n_truncated": n_trunc,
        "smallest_resolved_w": float(ws_clean.min()),
        "elapsed_s": dt,
    }
    with open("results/gaussian_family_positivity.json", "w") as f:
        json.dump(out, f, indent=2)
    print("Results: results/gaussian_family_positivity.json")


if __name__ == "__main__":
    main()
