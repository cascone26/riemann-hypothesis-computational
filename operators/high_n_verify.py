#!/usr/bin/env python3
"""
high_n_verify.py — N=512 verification for ambiguous cases + γ_eff analysis.

Key findings from full_prime_scan and discriminate_g6:
  γ₃: two stable j values — j=15(43) and j=22(251)
  γ₅: three stable j values — j=18(31), j=24(97), j=27(173)
  γ₆: two stable j values — j=29(127) and j=30(151), both stable through N=448
  γ₄: 491 drops to -9.874 at N=448 (oscillation dip?), 487 holds at -10.450

Discriminator: γ_eff = 2π×j/log(lsq) should equal γ_k for the true resonance.
Smaller |γ_eff - γ_k| = better match = true resonance.

Plan:
1. Compute γ_eff for all stable candidates; rank by error
2. Test N=512 for γ₃ (j=15/43 vs j=22/251), γ₅ (j=18/31 vs j=24/97 vs j=27/173)
3. Test N=512 for γ₆ (j=29/127 vs j=30/151) to see if 127 finally diverges
4. Test N=512 for γ₄ (491 vs 487) to confirm 491 oscillation
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
    return np.where(sv)[0]


primes5000 = make_sieve(5000)


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
    p_in = primes5000[primes5000 <= lsq]
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


ZEROS = {
    'g1': 14.134725141734693, 'g2': 21.022039638771555, 'g3': 25.010857580145688,
    'g4': 30.424876125859513, 'g5': 32.935061587739189, 'g6': 37.586178158651,
    'g7': 40.918719012680879, 'g8': 43.327073280914999, 'g9': 48.005150881167,
    'g10': 49.773832477672302
}


def gamma_eff(j, lsq):
    return TWO_PI * j / math.log(lsq)


def best_zero_match(j, lsq):
    ge = gamma_eff(j, lsq)
    best = min(ZEROS.items(), key=lambda kv: abs(kv[1] - ge))
    return best[0], best[1], abs(best[1] - ge)


# ===================================================================
# γ_eff ANALYSIS — compute best-zero match for all stable candidates
# ===================================================================
print("=" * 70)
print("γ_eff DISCRIMINATOR ANALYSIS")
print("True resonance = smallest |γ_eff - γ| among stable j values for each zero")
print("=" * 70)

candidates = [
    # (zero_name, j, lsq, label)
    ("γ₃", 15, 43, "near-miss?"),
    ("γ₃", 22, 251, "true?"),
    ("γ₅", 18, 31, "true? (best fit)"),
    ("γ₅", 24, 97, "competitor"),
    ("γ₅", 27, 173, "competitor"),
    ("γ₆", 29, 127, "near-miss?"),
    ("γ₆", 30, 151, "true?"),
    ("γ₇", 37, 293, "CONFIRMED UNIQUE"),
    ("γ₄", 30, 487, "competitor"),
    ("γ₄", 30, 491, "true?"),
]

print(f"\n{'Zero':6} {'j':4} {'lsq':6} {'γ_eff':8} {'|error|':8} {'best_zero':10} {'label'}")
print("-" * 70)
for (zname, j, lsq, label) in candidates:
    ge = gamma_eff(j, lsq)
    bname, bval, berr = best_zero_match(j, lsq)
    print(f"{zname:6} {j:4} {lsq:6} {ge:8.4f} {berr:8.4f} {bname:10} {label}")

# ===================================================================
# N=512 TESTS
# Priority: confirm which of the competing stable values diverges
# ===================================================================
print("\n\n" + "=" * 70)
print("N=512 VERIFICATION")
print("=" * 70)

test_cases = [
    # (label, lsq, known_trajectory)
    ("γ₃ j=15 lsq=43 (near-miss?)", 43, "N64=-11.72, N128=-11.78, N192=-11.22, N256=-10.65"),
    ("γ₃ j=22 lsq=251 (true?)", 251, "N64=-12.06, N128=-11.10, N192=-10.89, N256=-10.73"),
    ("γ₅ j=18 lsq=31 (true? best fit)", 31, "N64=-11.74, N128=-11.16, N192=-10.78, N256=-11.23"),
    ("γ₅ j=24 lsq=97 (competitor)", 97, "N64=-11.49, N128=-11.92, N192=-10.88, N256=-10.97"),
    ("γ₅ j=27 lsq=173 (competitor)", 173, "N64=-10.96, N128=-11.07, N192=-10.53, N256=-10.58"),
    ("γ₆ j=29 lsq=127 (near-miss?)", 127, "N256=-10.54, N320=-10.616, N448=-10.234"),
    ("γ₆ j=30 lsq=151 (true?)", 151, "N256=-11.54, N320=-10.568, N448=-10.261"),
    ("γ₄ j=30 lsq=487 (competitor)", 487, "N256=stable, N448=-10.450"),
    ("γ₄ j=30 lsq=491 (true?)", 491, "N256=stable, N448=-9.874 OSCILLATION?"),
]

print()
for label, lsq, prior in test_cases:
    print(f"\n--- {label} ---")
    print(f"  Prior: {prior}")
    e = compute_ev0(lsq, 512)
    log_e = math.log10(abs(e)) if abs(e) > 0 else -99
    ge_j = TWO_PI * round(sum(ZEROS[k] * math.log(lsq) / TWO_PI for k in ZEROS) / len(ZEROS)) / math.log(lsq)
    status = "STABLE" if log_e < -10 else "DIVERGING"
    print(f"  N=512: {log_e:.4f} → {status}", flush=True)

# ===================================================================
# SYNTHESIS
# ===================================================================
print("\n\n" + "=" * 70)
print("SYNTHESIS — which candidates survive N=512?")
print("=" * 70)
survivors_512 = {}
for label, lsq, prior in test_cases:
    e = compute_ev0(lsq, 512)
    log_e = math.log10(abs(e)) if abs(e) > 0 else -99
    survivors_512[label] = log_e
    status = "STABLE" if log_e < -10 else f"DIV ({log_e:.3f})"
    print(f"  {label[:40]:40} N512={log_e:7.3f} {status}")

print("\n\nAll done.")
