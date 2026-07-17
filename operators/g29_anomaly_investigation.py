#!/usr/bin/env python3
"""
g29_anomaly_investigation.py

The γ₂₉ precision-vs-stability paradox:
  lsq=15731: γ_eff err=0.000007 (BEST of all 30 zeros)  → ev₀(N256)=-8.24 (DIVERGES)
  lsq=19037: γ_eff err=0.000271 (40× worse)             → ev₀(N256)=-11.16 (STABLE)

Questions:
  1. Is lsq=15731 oscillating (like γ₂₀'s 4-cycle)? → test N=320..512
  2. What structural features differ between 15731 and 19037?
     - Prime count below each
     - Prime gap / next prime
     - Do any OTHER zeros resonate better with 15731?
  3. Systematic: among all top-4 candidates for zeros 21-30, do the "divergers"
     share a structural feature the stable ones don't?
"""
import math, numpy as np, scipy.linalg
from scipy.special import digamma

TWO_PI = 2 * math.pi
LOGPI = math.log(math.pi)
G_PER_N = 900

ZEROS_21_30 = [
    (21, 79.337375020249),
    (22, 82.910380854787),
    (23, 84.735492981329),
    (24, 87.425274613126),
    (25, 88.809111207634),
    (26, 92.491899270534),
    (27, 94.651344040819),
    (28, 95.870634228246),
    (29, 98.831194218008),
    (30, 101.317851006167),
]


def make_sieve(MAX_LSQ):
    sv = np.ones(MAX_LSQ + 1, dtype=bool)
    sv[0] = sv[1] = False
    for i in range(2, int(MAX_LSQ ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    return sv, np.where(sv)[0]


sv25000, primes25000 = make_sieve(25000)


def compute_ev0(lsq, N):
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
    p_in = primes25000[primes25000 <= lsq]
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
        return 0.0
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
    return float(scipy.linalg.eigvalsh(tau, subset_by_index=[0, 0])[0])


def gamma_eff_error(p, gamma):
    log_p = math.log(p)
    j = round(gamma * log_p / TWO_PI)
    ge = TWO_PI * j / log_p
    return abs(ge - gamma), j, ge


def best_resonance(lsq, zeros):
    """For a given prime lsq, what is the BEST γ_eff match among all 30 zeros?"""
    log_p = math.log(lsq)
    best_err = 1e9
    best_zero = None
    for k, gamma in zeros:
        j = round(gamma * log_p / TWO_PI)
        ge = TWO_PI * j / log_p
        err = abs(ge - gamma)
        if err < best_err:
            best_err = err
            best_zero = (k, gamma, j, ge)
    return best_err, best_zero


def prime_count(lsq):
    """Number of primes ≤ lsq."""
    return int(np.sum(primes25000 <= lsq))


def prime_gap(lsq):
    """Distance to next prime above lsq and gap size."""
    idx = np.searchsorted(primes25000, lsq)
    if idx < len(primes25000) and primes25000[idx] == lsq:
        # lsq is prime; find next prime
        if idx + 1 < len(primes25000):
            next_p = int(primes25000[idx + 1])
            prev_p = int(primes25000[idx - 1]) if idx > 0 else None
            gap_before = lsq - prev_p if prev_p else None
            gap_after = next_p - lsq
            return prev_p, lsq, next_p, gap_before, gap_after
    return None


# ============================================================
# PART 1: Is lsq=15731 oscillating or truly diverging?
# Prior: N64=-11.08, N128=-11.46, N192=-11.20, N256=-8.24
# ============================================================
print("=" * 70)
print("PART 1: lsq=15731 for γ₂₉ — Oscillation or True Divergence?")
print("Prior: N64=-11.08, N128=-11.46, N192=-11.20, N256=-8.24 (rose to -8!)")
print("Testing N=320..512 to see if it recovers below -10")
print("=" * 70)

GAMMA29 = 98.831194218008
P_ANOMALY = 15731

print(f"\nγ₂₉ = {GAMMA29}")
err_15731, j_15731, ge_15731 = gamma_eff_error(P_ANOMALY, GAMMA29)
print(f"lsq=15731: j={j_15731}, γ_eff={ge_15731:.6f}, err={err_15731:.6f}")
err_19037, j_19037, ge_19037 = gamma_eff_error(19037, GAMMA29)
print(f"lsq=19037: j={j_19037}, γ_eff={ge_19037:.6f}, err={err_19037:.6f}")

print("\nExtended N sweep for lsq=15731 (cap N=384; N=448/512 too large with Ollama at 10GB):")
for N in [64, 128, 192, 256, 320, 384]:
    e = compute_ev0(P_ANOMALY, N)
    log_e = math.log10(abs(e)) if abs(e) > 0 else -99
    src = "(prior)" if N <= 256 else "(new)"
    status = "stable" if log_e < -10 else "ABOVE"
    print(f"  N={N:4d}: log10|ev0|={log_e:>8.4f}  {status}  {src}", flush=True)

# ============================================================
# PART 2: Structural features — what distinguishes stable vs diverging?
# ============================================================
print("\n\n" + "=" * 70)
print("PART 2: Structural Feature Comparison — Stable vs Diverging Primes")
print("=" * 70)

# Known confirmed (stable) vs known diverging (but best γ_eff error)
# From zeros 21-30 scan:
#   γ₂₅: stable=409 (err=0.000324), diverged=18661 (err=0.000306) — smaller prime wins
#   γ₂₆: stable=281 (err=0.000399), larger candidates likely diverged
#   γ₂₉: stable=19037 (err=0.000271), diverged=15731 (err=0.000007)

STUDY_PAIRS = [
    # (gamma_k, gamma_val, stable_lsq, diverging_lsq, label)
    (25, 88.809111207634, 409, 18661, "γ₂₅: small vs large"),
    (29, 98.831194218008, 19037, 15731, "γ₂₉: precision paradox"),
]

for k, gamma, stable_p, div_p, label in STUDY_PAIRS:
    print(f"\n  {label}")
    for p, ptype in [(stable_p, "STABLE"), (div_p, "DIVERGED")]:
        info = prime_gap(p)
        n_primes = prime_count(p)
        err_r, _, _ = gamma_eff_error(p, gamma)

        # Best resonance with ANY of the 30 zeros
        best_err_r, best_zero_r = best_resonance(p, ZEROS_21_30)

        if info:
            prev_p, this_p, next_p, gap_before, gap_after = info
            gap_str = f"gap_before={gap_before}, gap_after={gap_after}, next={next_p}"
        else:
            gap_str = "not prime or out of range"

        print(f"    lsq={p:6d} ({ptype}): err={err_r:.6f}, π(lsq)={n_primes}, {gap_str}")
        print(f"                     best_zero=γ{best_zero_r[0]} (err={best_err_r:.6f})")


# ============================================================
# PART 3: For each pair (stable, best-err-diverger), compute ev0
# at a fine N grid to map the TRAJECTORY
# ============================================================
print("\n\n" + "=" * 70)
print("PART 3: Full N Trajectory — γ₂₉ Stable vs Diverging")
print("=" * 70)

NVALS_FINE = [32, 64, 96, 128, 160, 192, 224, 256, 320, 384]

for p, ptype in [(P_ANOMALY, "15731 DIVERGES"), (19037, "19037 STABLE")]:
    err_r, j_r, ge_r = gamma_eff_error(p, GAMMA29)
    print(f"\n  lsq={p} ({ptype}), err={err_r:.6f}")
    print(f"  N:  ", end="")
    vals = []
    for N in NVALS_FINE:
        e = compute_ev0(p, N)
        log_e = math.log10(abs(e)) if abs(e) > 0 else -99
        vals.append((N, log_e))
        print(f"N{N}={log_e:.3f}  ", end="", flush=True)
    print()
    above_10 = [(N, v) for N, v in vals if v > -10]
    if above_10:
        print(f"  Above -10 at: {[(N, round(v, 2)) for N, v in above_10]}")
    else:
        print(f"  All N values below -10: CONFIRMED STABLE")


# ============================================================
# PART 4: Is the anomaly connected to "competing resonances"?
# If lsq=15731 also resonates well with OTHER zeros, it might be "confused."
# ============================================================
print("\n\n" + "=" * 70)
print("PART 4: Competing Resonances — Does lsq=15731 also match other zeros?")
print("=" * 70)

ALL_30_ZEROS = [
    (1, 14.134725), (2, 21.022040), (3, 25.010858), (4, 30.424876), (5, 32.935062),
    (6, 37.586178), (7, 40.918719), (8, 43.327073), (9, 48.005151), (10, 49.773832),
    (11, 52.970321), (12, 56.446248), (13, 59.347044), (14, 60.831779), (15, 65.112544),
    (16, 67.079811), (17, 69.546402), (18, 72.067158), (19, 75.704691), (20, 77.144840),
    (21, 79.337375), (22, 82.910381), (23, 84.735493), (24, 87.425275), (25, 88.809111),
    (26, 92.491899), (27, 94.651344), (28, 95.870634), (29, 98.831194), (30, 101.317851),
]

for test_p in [15731, 19037]:
    log_p = math.log(test_p)
    print(f"\n  lsq={test_p}: top 6 resonances among first 30 zeros:")
    matches = []
    for k, gamma in ALL_30_ZEROS:
        j = round(gamma * log_p / TWO_PI)
        ge = TWO_PI * j / log_p
        err = abs(ge - gamma)
        matches.append((err, k, gamma, j, ge))
    matches.sort()
    for err, k, gamma, j, ge in matches[:6]:
        print(f"    γ{k:2d} ({gamma:.4f}): j={j}, γ_eff={ge:.4f}, err={err:.6f}")


print("\n\nDone. Key questions answered:")
print("  1. Is lsq=15731 oscillating or diverging? (see Part 1+3 trajectory)")
print("  2. Do competing resonances explain instability? (see Part 4)")
print("  3. Are prime-structural features (gap, count) predictive? (see Part 2)")
