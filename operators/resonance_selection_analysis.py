#!/usr/bin/env python3
"""
resonance_selection_analysis.py

THE PRECISION-VS-STABILITY PARADOX — ROOT CAUSE ANALYSIS

Observed in zeros 21-40 (dominant pattern):
  The prime with the SMALLEST γ_eff error is often NOT the stable one.
  Examples:
    γ₃₄: err=0.000070 (p=5441) diverges; err=0.000559 (p=19997) STABLE***
    γ₃₃: err=0.000100 (p=4127) diverges; err=0.000137 (p=5867) STABLE***
    γ₃₆: err=0.000085 (p=5009) borderline; err=0.000609 (p=773) STABLE***

HYPOTHESIS: γ_eff error = LOCAL resonance (single frequency match).
            Stability = GLOBAL compatibility (entire prime spectrum aligns).

The P-matrix P[n,m] = (S_m-S_n)/(π(m-n)) encodes the cumulative prime
distribution below lsq. For stability, P must closely approximate arch+B
across ALL matrix entries, not just the fundamental frequency.

PROPOSED NEW CRITERION: ||P - (arch+B)||_F / ||arch+B||_F (relative Frobenius norm)
This measures how well the ENTIRE prime spectrum cancels the zero-spectral operator.

If this Frobenius ratio is a better predictor of stability than γ_eff error,
that would be a fundamental theoretical result.
"""
import math, numpy as np, scipy.linalg
from scipy.special import digamma

TWO_PI = 2 * math.pi
LOGPI = math.log(math.pi)
G_PER_N = 900


def make_sieve(MAX_LSQ):
    sv = np.ones(MAX_LSQ + 1, dtype=bool)
    sv[0] = sv[1] = False
    for i in range(2, int(MAX_LSQ ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    return sv, np.where(sv)[0]


sv20000, primes20000 = make_sieve(20000)


def compute_all(lsq, N):
    """Return ev₀, and also arch, B, P matrices for structural analysis."""
    G = G_PER_N * N
    L = 2 * math.log(math.sqrt(lsq))
    ns = np.arange(-N, N + 1)
    D = len(ns)
    R = TWO_PI * N / L + 50.0
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h)
    w[0] = w[-1] = 0.5 * h
    Om = np.real(digamma(0.25 + 0.5j * grid)) - LOGPI
    c = TWO_PI * ns[:, None] / L
    d = c - grid[None, :]
    sm = np.abs(d) < 1e-12
    E = np.where(sm, np.sqrt(L) + 0j,
                 (1 / np.sqrt(L)) * np.exp(0.5j * d * L) * 2 * np.sin(0.5 * d * L) / d)
    arch = (E * (w * Om) @ E.conj().T).real / TWO_PI
    sh = np.sinh(L / 4) ** 2
    ni = ns[:, None]
    nj = ns[None, :]
    B = 32 * L * sh * (L ** 2 - 16 * math.pi ** 2 * ni * nj) / (
        (L ** 2 + 16 * math.pi ** 2 * ni ** 2) * (L ** 2 + 16 * math.pi ** 2 * nj ** 2))
    p_in = primes20000[primes20000 <= lsq]
    ll = []
    ww = []
    for p in p_in:
        lp = math.log(p)
        pk = p
        while pk <= lsq:
            ll.append(math.log(pk))
            ww.append(lp / math.sqrt(pk))
            pk *= p
    if not ll:
        return None
    lkpp = np.array(ll)
    wts_p = np.array(ww)
    K = len(lkpp)
    S = np.zeros(D)
    Sd = np.zeros(D)
    for k0 in range(0, K, 10000):
        k1 = min(k0 + 10000, K)
        ph = TWO_PI * np.outer(ns, lkpp[k0:k1]) / L
        S += np.sin(ph) @ wts_p[k0:k1]
        Sd += np.cos(ph) @ (2 * wts_p[k0:k1] * (1 - lkpp[k0:k1] / L))
    dn = ns[:, None] - ns[None, :]
    Si = S[:, None]
    Sj = S[None, :]
    with np.errstate(divide='ignore', invalid='ignore'):
        P = np.where(dn != 0, (Sj - Si) / (math.pi * dn), 0.0)
    np.fill_diagonal(P, Sd)
    tau = (arch + B - P)
    tau = (tau + tau.T) / 2
    ev0 = float(scipy.linalg.eigvalsh(tau, subset_by_index=[0, 0])[0])

    # Frobenius metrics
    AB = arch + B
    frob_P = np.linalg.norm(P, 'fro')
    frob_AB = np.linalg.norm(AB, 'fro')
    frob_diff = np.linalg.norm(P - AB, 'fro')
    frob_ratio = frob_diff / frob_AB if frob_AB > 0 else 999

    return ev0, frob_P, frob_AB, frob_diff, frob_ratio


def gamma_eff_error(p, gamma):
    log_p = math.log(p)
    j = round(gamma * log_p / TWO_PI)
    ge = TWO_PI * j / log_p
    return abs(ge - gamma), j, ge


# Test cases: (gamma_k, gamma_val, [(lsq, known_verdict)])
TEST_CASES = [
    (33, 107.168611184, [
        (4127, "DIVERGES"),    # err=0.000100 (best), diverges
        (5867, "STABLE***"),   # err=0.000137, stable
        (6221, "div(neg)"),    # err=0.000376, late-converger
        (953,  "STABLE(N256)"), # err=0.000397, borderline
    ]),
    (34, 111.029535543, [
        (5441, "DIVERGES"),    # err=0.000070 (best), diverges
        (6823, "DIVERGES"),    # err=0.000239, diverges
        (1399, "div(neg)"),    # err=0.000530, late-converger
        (19997,"STABLE***"),   # err=0.000559 (worst), stable
    ]),
    (36, 114.320220915, [
        (5009, "STABLE(N256)"), # err=0.000085 (best), borderline
        (5591, "DIVERGES"),    # err=0.000093, diverges
        (11423,"DIVERGES"),    # err=0.000403, diverges
        (773,  "STABLE***"),   # err=0.000609 (worst), stable
    ]),
    (38, 118.790782866, [
        (19751,"DIVERGES"),    # err=0.000085 (best), diverges
        (9419, "STABLE***"),   # err=0.000232, stable
        (15161,"DIVERGES"),    # err=0.000253, diverges
        (2381, "STABLE***"),   # err=0.000368, stable
    ]),
]

N_ANALYSIS = 128  # Fast analysis N (use 128 for matrix inspection)

print("=" * 80)
print("RESONANCE SELECTION ANALYSIS")
print("Testing hypothesis: Frobenius norm ||P - (arch+B)||_F predicts stability")
print("better than γ_eff frequency error")
print("=" * 80)
print()
print(f"{'Zero':>5} {'lsq':>8} {'γ_eff_err':>12} {'frob_ratio':>12} {'ev0(N128)':>12} {'verdict':>20}")
print("-" * 80)

all_results = []

for k, gamma, pairs in TEST_CASES:
    print(f"\n  γ{k} = {gamma:.6f}")
    for lsq, verdict in pairs:
        err, j, ge = gamma_eff_error(lsq, gamma)
        result = compute_all(lsq, N_ANALYSIS)
        if result is None:
            print(f"  {k:>5} {lsq:>8} {err:>12.6f}  {'N/A':>12}  {'N/A':>12} {verdict:>20}")
            continue
        ev0, frob_P, frob_AB, frob_diff, frob_ratio = result
        log_ev0 = math.log10(abs(ev0)) if abs(ev0) > 0 else -99
        is_stable = "STABLE" in verdict
        all_results.append({
            'k': k, 'gamma': gamma, 'lsq': lsq, 'err': err,
            'frob_ratio': frob_ratio, 'ev0': log_ev0, 'verdict': verdict,
            'is_stable': is_stable
        })
        star = "<<STABLE" if is_stable else ""
        print(f"    {k:>5} {lsq:>8} {err:>12.6f} {frob_ratio:>12.6f} {log_ev0:>12.4f} {verdict:>20}  {star}", flush=True)

print()
print("=" * 80)
print("CORRELATION ANALYSIS")
print("=" * 80)

# Rank by γ_eff error vs stability
print("\nRanked by γ_eff error (smaller=better precision):")
for k in set(r['k'] for r in all_results):
    cases = [r for r in all_results if r['k'] == k]
    cases.sort(key=lambda r: r['err'])
    print(f"  γ{k}:")
    for rank, r in enumerate(cases):
        flag = "← STABLE" if r['is_stable'] else ""
        print(f"    rank {rank+1}: lsq={r['lsq']:6d}, err={r['err']:.6f}, frob={r['frob_ratio']:.6f}, ev0={r['ev0']:.3f}  {flag}")

print("\nRanked by Frobenius ratio (smaller=better cancellation):")
for k in set(r['k'] for r in all_results):
    cases = [r for r in all_results if r['k'] == k]
    cases.sort(key=lambda r: r['frob_ratio'])
    print(f"  γ{k}:")
    for rank, r in enumerate(cases):
        flag = "← STABLE" if r['is_stable'] else ""
        print(f"    rank {rank+1}: lsq={r['lsq']:6d}, frob={r['frob_ratio']:.6f}, err={r['err']:.6f}, ev0={r['ev0']:.3f}  {flag}")

print()
print("KEY QUESTION: For each zero, does the SMALLEST Frobenius ratio")
print("correspond to the STABLE prime?")
print("If YES → Frobenius ratio is the true resonance criterion, not γ_eff error.")
print("\nDone.")
