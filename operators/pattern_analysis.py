#!/usr/bin/env python3
"""
pattern_analysis.py

Deep analysis of all 50 confirmed (γ_k, p_k, j_k) pairs.
Goals:
  1. Find the mathematical law governing p_k for γ_k
  2. Analyze j_k distribution: is j_k ~ γ_k * log(p_k) / (2π) rational?
  3. Asymptotic: does p_k grow with k? Does log(p_k)/log(γ_k) converge?
  4. Twin prime / prime gap analysis for stable vs unstable candidates
  5. Check if j_k / k has a pattern (arithmetic progressions?)
"""
import math
import numpy as np

TWO_PI = 2 * math.pi

# All 50 confirmed zeros — (k, gamma_k, p_k, j_k, err)
CONFIRMED_50 = [
    # zeros 1-10
    (1,  14.134725, 31,    2,   0.000396),
    (2,  21.022040, 199,   13,  0.000198),
    (3,  25.010858, 127,   14,  0.000168),
    (4,  30.424876, 491,   27,  0.005000),
    (5,  32.935062, 31,    4,   0.000400),
    (6,  37.586178, 151,   22,  0.017000),
    (7,  40.918719, 293,   28,  0.009000),
    (8,  43.327073, 79,    18,  0.000028),
    (9,  48.005151, 461,   38,  0.000350),
    (10, 49.773832, 709,   44,  0.004000),
    # zeros 11-20
    (11, 52.970321, 103,   32,  0.000142),
    (12, 56.446248, 167,   38,  0.026000),
    (13, 59.347044, 619,   54,  0.000212),
    (14, 60.831779, 277,   46,  0.000033),
    (15, 65.112544, 229,   51,  0.000072),
    (16, 67.079811, 641,   66,  0.000154),
    (17, 69.546402, 1511,  91,  0.000212),
    (18, 72.067158, 1973,  100, 0.000054),
    (19, 75.704691, 2311,  111, 0.000039),
    (20, 77.144840, 733,   83,  0.000343),
    # zeros 21-30
    (21, 79.337375, 19531, 191, 0.000015),
    (22, 82.910381, 2243,  115, 0.000080),
    (23, 84.735493, 2351,  118, 0.000110),
    (24, 87.425275, 10627, 162, 0.000238),
    (25, 88.809111, 409,   73,  0.000324),
    (26, 92.491899, 281,   77,  0.000399),
    (27, 94.651344, 18211, 187, 0.000140),
    (28, 95.870634, 17597, 186, 0.000007),
    (29, 98.831194, 19037, 197, 0.000271),
    (30, 101.317851, 18253, 202, 0.000065),
    # zeros 31-40
    (31, 103.725538, 4019,  148, 0.000012),
    (32, 105.446623, 2053,  139, 0.000036),
    (33, 107.168611, 5867,  165, 0.000137),
    (34, 111.029536, 19997, 228, 0.000559),
    (35, 111.874659, 5393,  160, 0.000511),
    (36, 114.320221, 5009,  165, 0.000085),
    (37, 116.226680, 18749, 237, 0.000250),
    (38, 118.790783, 9419,  206, 0.000232),
    (39, 121.370125, 313,   96,  0.002814),
    (40, 122.946829, 18253, 249, 0.000700),
    # zeros 41-50
    (41, 124.256819, 14143, 189, 0.000283),
    (42, 127.516684, 2179,  156, 0.000578),
    (43, 129.578704, 2341,  160, 0.000664),
    (44, 131.087689, 2141,  160, 0.000723),
    (45, 133.497737, 17029, 207, 0.000559),
    (46, 134.756510, 373,   127, 0.001133),
    (47, 138.116042, 1103,  154, 0.000187),
    (48, 139.736209, 9631,  204, 0.000606),
    (49, 141.123470, 5393,  193, 0.000111),
    (50, 143.111846, 2269,  176, 0.000241),
]


def make_sieve(MAX):
    sv = np.ones(MAX + 1, dtype=bool)
    sv[0] = sv[1] = False
    for i in range(2, int(MAX**0.5) + 1):
        if sv[i]:
            sv[i*i::i] = False
    return np.where(sv)[0]

primes = make_sieve(25000)


def prime_gap(p):
    idx = np.searchsorted(primes, p)
    if idx < len(primes) and primes[idx] == p:
        prev_p = int(primes[idx-1]) if idx > 0 else None
        next_p = int(primes[idx+1]) if idx+1 < len(primes) else None
        gb = p - prev_p if prev_p else None
        ga = next_p - p if next_p else None
        return gb, ga
    return None, None


print("="*80)
print("Pattern Analysis — 50 Confirmed (γ_k, p_k, j_k) Pairs")
print("="*80)
print()

# ============================================================
# 1. Basic statistics
# ============================================================
ks = np.array([x[0] for x in CONFIRMED_50])
gammas = np.array([x[1] for x in CONFIRMED_50])
ps = np.array([x[2] for x in CONFIRMED_50])
js = np.array([x[3] for x in CONFIRMED_50])
errs = np.array([x[4] for x in CONFIRMED_50])

print("1. BASIC STATISTICS")
print(f"   p_k range: [{ps.min()}, {ps.max()}]   median={np.median(ps):.0f}   mean={ps.mean():.0f}")
print(f"   j_k range: [{js.min()}, {js.max()}]   median={np.median(js):.0f}")
print(f"   err range: [{errs.min():.6f}, {errs.max():.6f}]")
print()

# ============================================================
# 2. Frequency ratio j/k — does it converge?
# ============================================================
print("2. HARMONIC NUMBER j_k/k:")
j_over_k = js / ks
print(f"   j_k/k range: [{j_over_k.min():.2f}, {j_over_k.max():.2f}]")
print(f"   j_k/k for k=1..10:  {[f'{v:.2f}' for v in j_over_k[:10]]}")
print(f"   j_k/k for k=41..50: {[f'{v:.2f}' for v in j_over_k[40:]]}")
print()

# ============================================================
# 3. Effective frequency: γ_k / j_k ~ 2π / log(p_k)
# ============================================================
print("3. γ_EFF RESIDUAL ANALYSIS (γ_k - 2π*j_k/log(p_k)):")
residuals = []
for k, gamma, p, j, err in CONFIRMED_50:
    log_p = math.log(p)
    ge = TWO_PI * j / log_p
    res = gamma - ge
    residuals.append(res)

residuals = np.array(residuals)
print(f"   Residual (γ_k - γ_eff): mean={residuals.mean():.6f}, std={residuals.std():.6f}")
print(f"   |residual| > 0.001: {sum(1 for r in residuals if abs(r) > 0.001)} cases")
print(f"   |residual| > 0.01:  {sum(1 for r in residuals if abs(r) > 0.01)} cases")
print()

# ============================================================
# 4. Log(p_k) / log(k) — does p_k ~ k^alpha?
# ============================================================
print("4. PRIME GROWTH: log(p_k)/log(k)  — does p_k ~ k^alpha?")
log_p = np.log(ps)
log_k = np.log(ks)
ratios = log_p / log_k
print(f"   log(p_k)/log(k): mean={ratios.mean():.3f}, std={ratios.std():.3f}")
print(f"   First 10: {[f'{v:.2f}' for v in ratios[:10]]}")
print(f"   Last 10:  {[f'{v:.2f}' for v in ratios[40:]]}")
print()

# ============================================================
# 5. j_k vs gamma_k — does j ~ gamma * C for some C?
# ============================================================
print("5. HARMONIC SCALING: j_k / γ_k  — if p_k is the same for many zeros...")
j_over_gamma = js / gammas
print(f"   j_k/γ_k: mean={j_over_gamma.mean():.4f}, std={j_over_gamma.std():.4f}")
print(f"   This equals log(p_k)/(2π) on average = {np.log(ps).mean()/TWO_PI:.4f}")
print()

# ============================================================
# 6. Which primes are REUSED across different zeros?
# ============================================================
from collections import Counter
prime_counts = Counter(ps)
reused = {p: c for p, c in prime_counts.items() if c > 1}
print("6. PRIME REUSE — same p_k confirmed for multiple zeros:")
if reused:
    for p, c in sorted(reused.items()):
        zeros_using_p = [(k, gamma) for k, gamma, pk, j, err in CONFIRMED_50 if pk == p]
        print(f"   p={p}: used for γ_{zeros_using_p[0][0]}={zeros_using_p[0][1]:.4f}, γ_{zeros_using_p[1][0]}={zeros_using_p[1][1]:.4f}")
        if c > 2:
            for kz, gz in zeros_using_p[2:]:
                print(f"          also γ_{kz}={gz:.4f}")
else:
    print("   No prime reuse — all 50 zeros use distinct primes")
print()

# ============================================================
# 7. j_k repetitions — same harmonic number for different zeros?
# ============================================================
j_counts = Counter(js)
reused_j = {j: c for j, c in j_counts.items() if c > 1}
print("7. HARMONIC REPETITIONS — same j_k for multiple zeros:")
if reused_j:
    for j, c in sorted(reused_j.items()):
        zeros_using_j = [(k, gamma, p) for k, gamma, p, jj, err in CONFIRMED_50 if jj == j]
        print(f"   j={j}: ", end="")
        print(", ".join(f"γ_{kz}={gz:.4f}(p={pz})" for kz, gz, pz in zeros_using_j))
else:
    print("   No j repetitions")
print()

# ============================================================
# 8. Prime gap analysis — stable vs all confirmed
# ============================================================
print("8. PRIME GAP ANALYSIS for confirmed stable primes:")
gap_2 = gap_4 = gap_6 = gap_other = 0
small_p_count = 0
for k, gamma, p, j, err in CONFIRMED_50:
    gb, ga = prime_gap(p)
    if ga == 2: gap_2 += 1
    elif ga == 4: gap_4 += 1
    elif ga == 6: gap_6 += 1
    else: gap_other += 1
    if p < 1000:
        small_p_count += 1
print(f"   gap_after=2 (twin):   {gap_2}/50 ({100*gap_2//50}%)")
print(f"   gap_after=4 (cousin): {gap_4}/50 ({100*gap_4//50}%)")
print(f"   gap_after=6 (sexy):   {gap_6}/50 ({100*gap_6//50}%)")
print(f"   gap_after>=8:         {gap_other}/50 ({100*gap_other//50}%)")
print(f"   Small primes p<1000:  {small_p_count}/50")
print()

# ============================================================
# 9. Key ratio: γ_k * log(p_k) = 2π * j_k → check
# ============================================================
print("9. FUNDAMENTAL RELATION CHECK: γ_k * log(p_k) ≈ 2π * j_k")
print("   (should hold by construction — checking precision and distribution)")
diffs = []
for k, gamma, p, j, err in CONFIRMED_50:
    lhs = gamma * math.log(p)
    rhs = TWO_PI * j
    diffs.append(lhs - rhs)
diffs = np.array(diffs)
print(f"   Difference γ*log(p) - 2π*j: mean={diffs.mean():.4f}, std={diffs.std():.4f}")
print(f"   Range: [{diffs.min():.4f}, {diffs.max():.4f}]")
print()

# ============================================================
# 10. Is there a Diophantine approximation pattern?
# j_k = round(γ_k * log(p_k) / 2π) — best rational approx to γ_k / (2π/log(p_k))
# ============================================================
print("10. DIOPHANTINE APPROXIMATION: quality of γ_k * log(p_k) / (2π) ≈ j_k")
print("    (convergents of continued fraction expansion)")
print()

# Check which confirmed primes give exceptionally small err
print("    Top 10 smallest γ_eff errors:")
sorted_by_err = sorted(CONFIRMED_50, key=lambda x: x[4])
for k, gamma, p, j, err in sorted_by_err[:10]:
    print(f"    γ_{k:2d}={gamma:.6f}: p={p:6d}, j={j:4d}, err={err:.7f}")

print()
print("    Top 10 largest γ_eff errors:")
for k, gamma, p, j, err in sorted_by_err[-10:]:
    print(f"    γ_{k:2d}={gamma:.6f}: p={p:6d}, j={j:4d}, err={err:.7f}")
print()

# ============================================================
# 11. Asymptotic: 2π*j_k / γ_k = log(p_k) → p_k ~ e^{γ_k * 2π/j_k}?
# But j_k varies. What IS j_k/γ_k?
# ============================================================
print("11. ASYMPTOTIC BEHAVIOR: log(p_k) vs γ_k")
log_p_vals = np.log(ps)
print(f"    Correlation(log(p_k), γ_k) = {np.corrcoef(log_p_vals, gammas)[0,1]:.4f}")
# Linear fit log(p_k) = a * gamma_k + b
A = np.vstack([gammas, np.ones(len(gammas))]).T
slope, intercept = np.linalg.lstsq(A, log_p_vals, rcond=None)[0]
print(f"    Linear fit: log(p_k) ≈ {slope:.4f} * γ_k + {intercept:.4f}")
print(f"    => p_k ~ exp({intercept:.4f}) * γ_k^{slope:.4f}  [if slope≈0, p_k is bounded]")
print(f"    => p_k ~ e^({intercept:.3f}) * e^({slope:.3f}*γ_k)")
print(f"    => p_k ~ {math.exp(intercept):.1f} * e^({slope:.4f}*γ_k)")
print()

# ============================================================
# 12. j_k sequence — is it approximately arithmetic?
# ============================================================
print("12. j_k SEQUENCE ANALYSIS:")
print(f"    j_k: {list(js[:20])}")
print(f"    j_k: {list(js[20:40])}")
print(f"    j_k: {list(js[40:])}")
dj = np.diff(js)
print(f"    Differences Δj_k: mean={dj.mean():.2f}, std={dj.std():.2f}, range=[{dj.min()}, {dj.max()}]")
print(f"    Δj_k sequence (first 20): {list(dj[:20])}")
print()

print("="*80)
print("CONCLUSIONS:")
# Prime reuse
if reused:
    print(f"  PRIME REUSE: {len(reused)} primes serve multiple zeros ({list(reused.items())})")
else:
    print("  PRIME REUSE: None — each zero gets a unique prime (50 distinct primes)")

# j distribution
print(f"  HARMONIC RANGE: j_k ∈ [{js.min()}, {js.max()}] (ratio max/min = {js.max()/js.min():.1f})")

# Growth
print(f"  PRIME GROWTH: slope={slope:.4f} → p_k {'grows exponentially with γ' if slope > 0.05 else 'does NOT grow strongly with γ'}")
print(f"  CORRELATION: log(p_k)~γ_k correlation = {np.corrcoef(log_p_vals, gammas)[0,1]:.4f}")
print("="*80)
