#!/usr/bin/env python3
"""
twin_prime_boundary_test.py

Validates the Twin Prime Boundary Hypothesis:
  "If the best-γ_eff-error candidate for zero γ_k diverges,
   does it have gap_after ≤ 4 (twin/cousin prime at boundary)?"

Tests this across all 40 confirmed zeros by:
  1. Finding the best-err candidate for each zero
  2. Checking its gap_after
  3. Comparing to whether it was the confirmed stable prime or the diverging one

Data: all 40 confirmed stable primes (from scan results) and best-err alternatives.
"""
import math
import numpy as np

TWO_PI = 2 * math.pi


def make_sieve(MAX_LSQ):
    sv = np.ones(MAX_LSQ + 1, dtype=bool)
    sv[0] = sv[1] = False
    for i in range(2, int(MAX_LSQ ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    return sv, np.where(sv)[0]


sv, primes = make_sieve(30000)


def prime_gap(p, primes):
    """Return (gap_before, gap_after) for prime p."""
    idx = np.searchsorted(primes, p)
    if idx < len(primes) and primes[idx] == p:
        prev_p = int(primes[idx - 1]) if idx > 0 else None
        next_p = int(primes[idx + 1]) if idx + 1 < len(primes) else None
        gb = p - prev_p if prev_p else None
        ga = next_p - p if next_p else None
        return gb, ga
    return None, None


def best_candidate(gamma, p_min=50, p_max=20000):
    """Find best γ_eff error candidate for given zero."""
    best = None
    for p in primes:
        if p < p_min or p > p_max:
            continue
        log_p = math.log(p)
        j = round(gamma * log_p / TWO_PI)
        if j < 1:
            continue
        ge = TWO_PI * j / log_p
        err = abs(ge - gamma)
        if best is None or err < best[0]:
            best = (err, j, p, ge)
    return best


def top_candidates(gamma, p_min=50, p_max=20000, top_n=4):
    """Find top-N candidates by γ_eff error."""
    cands = []
    for p in primes:
        if p < p_min or p > p_max:
            continue
        log_p = math.log(p)
        j = round(gamma * log_p / TWO_PI)
        if j < 1:
            continue
        ge = TWO_PI * j / log_p
        err = abs(ge - gamma)
        cands.append((err, j, p, ge))
    cands.sort()
    return cands[:top_n]


# Confirmed stable primes from scan results (zeros 1-40)
# Format: (zero_k, gamma, stable_lsq, stable_err, paradox=True if best-err diverges)
CONFIRMED = [
    # zeros 1-10 (from paper)
    (1,  14.134725, 31,    0.000396, False),  # γ₁: small j=2, best-err wins
    (2,  21.022040, 199,   0.000198, True),   # γ₂: 199 confirmed, best might differ
    (3,  25.010858, 127,   0.000168, False),
    (4,  30.424876, 491,   0.005000, True),   # γ₄: best-err lsq might differ
    (5,  32.935062, 31,    0.000400, False),  # γ₅: lsq=31, j=4
    (6,  37.586178, 151,   0.017000, False),
    (7,  40.918719, 293,   0.009000, False),
    (8,  43.327073, 79,    0.000028, False),
    (9,  48.005151, 461,   0.000350, False),
    (10, 49.773832, 709,   0.004000, False),
    # zeros 11-20
    (11, 52.970321, 103,   0.000142, False),
    (12, 56.446248, 167,   0.026000, True),   # γ₁₂: confirmed 167, best-err (4721) diverges
    (13, 59.347044, 619,   0.000212, False),
    (14, 60.831779, 277,   0.000033, False),
    (15, 65.112544, 229,   0.000072, False),
    (16, 67.079811, 641,   0.000154, False),
    (17, 69.546402, 1511,  0.000212, False),
    (18, 72.067158, 1973,  0.000054, False),
    (19, 75.704691, 2311,  0.000039, False),
    (20, 77.144840, 733,   0.000343, False),  # γ₂₀: dip at N128-192 but recovers
    # zeros 21-30
    (21, 79.337375, 19531, 0.000015, False),
    (22, 82.910381, 2243,  0.000080, False),
    (23, 84.735493, 2351,  0.000110, False),
    (24, 87.425275, 10627, 0.000238, False),
    (25, 88.809111, 409,   0.000324, True),   # γ₂₅: 409 stable; larger candidates diverge
    (26, 92.491899, 281,   0.000399, True),   # γ₂₆: small prime stable; larger diverge
    (27, 94.651344, 18211, 0.000140, True),   # Need to check
    (28, 95.870634, 17597, 0.000007, False),  # γ₂₈: best-err is stable here!
    (29, 98.831194, 19037, 0.000271, True),   # γ₂₉: 15731 (err=7e-6) DIVERGES
    (30, 101.317851, 18253, 0.000065, False),
    # zeros 31-40
    (31, 103.725538, 4019,  0.000012, False),
    (32, 105.446623, 2053,  0.000036, True),  # γ₃₂: late converger; best-err may differ
    (33, 107.168611, 5867,  0.000137, True),  # γ₃₃: 4127 (best-err) diverges
    (34, 111.029536, 19997, 0.000559, True),  # γ₃₄: top-3 diverge
    (35, 111.874659, 5393,  0.000511, True),  # γ₃₅: top-2 diverge
    (36, 114.320221, 5009,  0.000085, True),  # γ₃₆: best (11423) diverges
    (37, 116.226680, 18749, 0.000250, False),
    (38, 118.790783, 9419,  0.000232, True),  # γ₃₈: best (19751) diverges
    (39, 121.370125, 313,   0.002814, True),  # γ₃₉: all err<0.001 diverge!
    (40, 122.946829, 18253, 0.000700, False),
]


print("=" * 80)
print("Twin Prime Boundary Hypothesis Validation")
print("Question: Do diverging best-err candidates tend to have small gap_after?")
print("=" * 80)
print()

paradox_cases = []
no_paradox_cases = []

for k, gamma, stable_lsq, stable_err, paradox in CONFIRMED:
    # Find true best-err candidate in [50, 20000]
    top4 = top_candidates(gamma, p_min=50, p_max=20000, top_n=4)
    best_err, best_j, best_p, best_ge = top4[0]

    gb_stable, ga_stable = prime_gap(stable_lsq, primes)
    gb_best, ga_best = prime_gap(best_p, primes)

    is_twin_stable = ga_stable is not None and ga_stable <= 4
    is_twin_best = ga_best is not None and ga_best <= 4

    if paradox and best_p != stable_lsq:
        # True paradox case: best-err diverges, stable is different
        paradox_cases.append({
            'k': k, 'gamma': gamma,
            'stable_lsq': stable_lsq, 'stable_err': stable_err,
            'ga_stable': ga_stable, 'is_twin_stable': is_twin_stable,
            'best_p': best_p, 'best_err': best_err,
            'ga_best': ga_best, 'is_twin_best': is_twin_best,
        })
    else:
        no_paradox_cases.append({
            'k': k, 'gamma': gamma,
            'stable_lsq': stable_lsq, 'stable_err': stable_err,
            'ga_stable': ga_stable, 'is_twin_stable': is_twin_stable,
            'best_p': best_p, 'best_err': best_err,
            'ga_best': ga_best, 'is_twin_best': is_twin_best,
        })

print(f"PARADOX CASES ({len(paradox_cases)} zeros — best-err candidate diverges):")
print(f"{'Zero':<6} {'γ_k':<10} {'best_p(div)':<14} {'best_err':<10} {'gap_after':<10} {'twin?':<6} {'stable_p':<10} {'gap_after_stable':<17} {'twin?'}")
print("-" * 95)
twin_divergers = 0
for c in paradox_cases:
    ga_best_str = str(c['ga_best']) if c['ga_best'] else "?"
    ga_stable_str = str(c['ga_stable']) if c['ga_stable'] else "?"
    twin_best = "TWIN" if c['is_twin_best'] else "-"
    twin_stable = "twin" if c['is_twin_stable'] else "-"
    if c['is_twin_best']:
        twin_divergers += 1
    print(f"γ{c['k']:<5} {c['gamma']:<10.6f} {c['best_p']:<14} {c['best_err']:<10.6f} {ga_best_str:<10} {twin_best:<6} {c['stable_lsq']:<10} {ga_stable_str:<17} {twin_stable}")

print()
print(f"  Of {len(paradox_cases)} paradox cases, {twin_divergers} ({100*twin_divergers//max(1,len(paradox_cases))}%) have gap_after ≤ 4 on the diverging candidate")

print()
print(f"NON-PARADOX CASES ({len(no_paradox_cases)} zeros — best-err candidate stable):")
print(f"{'Zero':<6} {'γ_k':<10} {'stable_p':<14} {'err':<10} {'gap_after':<10} {'twin?'}")
print("-" * 60)
twin_stable_count = 0
for c in no_paradox_cases:
    ga_str = str(c['ga_stable']) if c['ga_stable'] else "?"
    twin_str = "TWIN" if c['is_twin_stable'] else "-"
    if c['is_twin_stable']:
        twin_stable_count += 1
    print(f"γ{c['k']:<5} {c['gamma']:<10.6f} {c['stable_lsq']:<14} {c['stable_err']:<10.6f} {ga_str:<10} {twin_str}")

print()
print(f"  Of {len(no_paradox_cases)} non-paradox cases, {twin_stable_count} ({100*twin_stable_count//max(1,len(no_paradox_cases))}%) have gap_after ≤ 4")

print()
print("=" * 80)
print("HYPOTHESIS PREDICTION: High rate of gap_after ≤ 4 for DIVERGING best-err candidates,")
print("                       low rate for stable confirmed primes.")
print("=" * 80)

# Additional: check the specific γ₂₉ case in detail
print()
print("γ₂₉ detail:")
gb15731, ga15731 = prime_gap(15731, primes)
gb19037, ga19037 = prime_gap(19037, primes)
print(f"  lsq=15731 (DIVERGES): gap_before={gb15731}, gap_after={ga15731} → twin={ga15731 <= 4}")
print(f"  lsq=19037 (STABLE):   gap_before={gb19037}, gap_after={ga19037} → twin={ga19037 <= 4}")
