#!/usr/bin/env python3
"""
fake_zero_test.py

THE KEY TEST for RH connection:
  For real Riemann zeros γ_k, the CvS operator finds a prime lsq with ev0 → -∞ (PSD)
  If this is connected to RH, it should FAIL for fake γ values (not true zeros).

We test γ_fake ∉ {true Riemann zeros} and check if any prime produces STABLE.
If NO fake γ produces a stable lsq → operator genuinely detects critical-line zeros.

Fake values chosen to span different ranges, including values near but distinct from true zeros.
"""
import math, numpy as np, scipy.linalg, sys
from scipy.special import digamma

TWO_PI = 2 * math.pi
LOGPI = math.log(math.pi)
G_PER_N = 900

# True zeros (reference) — first 10
TRUE_ZEROS_10 = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
                 37.586178, 40.918719, 43.327073, 48.005151, 49.773832]

# Fake γ values — NOT Riemann zeros, testing operator discrimination
# Chosen to be: between zeros, near zeros but offset, round numbers in zero range
FAKE_ZEROS = [
    (1, 13.5,    "below γ₁"),
    (2, 14.0,    "just below γ₁=14.135"),
    (3, 14.5,    "just above γ₁"),
    (4, 15.0,    "between γ₁ and γ₂"),
    (5, 18.0,    "between γ₁ and γ₂"),
    (6, 20.0,    "just below γ₂=21.022"),
    (7, 21.5,    "just above γ₂"),
    (8, 25.5,    "just above γ₃=25.011"),
    (9, 30.0,    "just below γ₄=30.425"),
    (10, 40.0,   "between γ₆ and γ₇"),
    (11, 50.5,   "just above γ₁₀=49.774"),
    (12, 98.0,   "near γ₂₉=98.831 but different"),
    (13, 99.5,   "between γ₂₉ and γ₃₀"),
    (14, 143.0,  "near γ₅₀=143.112 but different"),
]


def make_sieve(MAX_LSQ):
    sv = np.ones(MAX_LSQ + 1, dtype=bool)
    sv[0] = sv[1] = False
    for i in range(2, int(MAX_LSQ ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    return sv, np.where(sv)[0]


sv20000, primes20000 = make_sieve(20000)


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


def prime_candidates(gamma, p_min=50, p_max=20000, top_n=12):
    cands = []
    for p in primes20000:
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


def test_fake_zero(fake_id, gamma, description):
    print(f"\n{'='*70}")
    print(f"FAKE γ#{fake_id} = {gamma:.4f}  [{description}]")
    print(f"{'='*70}")
    sys.stdout.flush()

    cands = prime_candidates(gamma, top_n=12)
    best_n256 = None
    any_stable = False

    for err, j, p, ge in cands[:12]:
        n_max = 256 if p >= 5000 else 320
        nvals = [64, 128, 192, n_max]
        row = {}
        for N in nvals:
            e = compute_ev0(p, N)
            row[N] = math.log10(abs(e)) if abs(e) > 0 else -99
        n_last = nvals[-1]
        delta = row[n_last] - row[64]
        stable = (row[n_last] < -10) and (delta < 2.5)
        all_stable = all(row[N] < -10 for N in nvals)
        vals = "  ".join(f"N{N}={row[N]:>7.4f}" for N in nvals)
        if all_stable:
            verdict = "STABLE***"
        elif stable:
            verdict = f"STABLE(1side,Δ={delta:+.2f})"
        else:
            verdict = f"div(Δ={delta:+.2f})"
        print(f"  p={p:6d} j={j:4d} err={err:.6f}: {vals}  → {verdict}", flush=True)
        if stable:
            any_stable = True
        if best_n256 is None or row[n_last] < best_n256[1]:
            best_n256 = (p, row[n_last], err, verdict)

    print(f"\n  Best ev0 at N_max: p={best_n256[0]}, log10|ev0|={best_n256[1]:.4f}, err={best_n256[2]:.6f}")
    if any_stable:
        print(f"  *** STABLE FOUND for fake γ={gamma} — operator does NOT discriminate! ***")
    else:
        print(f"  No stable lsq found — operator correctly rejects fake zero (best ev0={best_n256[1]:.4f})")
    return any_stable, best_n256[1]


print("="*70)
print("FAKE ZERO TEST — RH Connection Verification")
print("Do true Riemann zeros have STABLE lsq that fake zeros do NOT?")
print("="*70)
print()

results = []
for fake_id, gamma, desc in FAKE_ZEROS:
    stable_found, best_ev0 = test_fake_zero(fake_id, gamma, desc)
    results.append((fake_id, gamma, desc, stable_found, best_ev0))
    sys.stdout.flush()

print("\n\n" + "="*70)
print("SUMMARY — Fake Zero Discrimination Test")
print("="*70)
n_discriminated = sum(1 for _, _, _, s, _ in results if not s)
n_failed = sum(1 for _, _, _, s, _ in results if s)
for fake_id, gamma, desc, stable_found, best_ev0 in results:
    status = "REJECTED (correct!)" if not stable_found else "HAS STABLE (operator fails!)"
    print(f"  γ_fake={gamma:<8.4f} [{desc:<35}]: {status}  best_ev0={best_ev0:.4f}")

print()
print(f"Operator correctly rejected: {n_discriminated}/{len(results)} fake zeros")
print(f"Operator failed to reject:   {n_failed}/{len(results)} fake zeros")
if n_failed == 0:
    print("RESULT: Operator PERFECTLY discriminates true vs fake zeros.")
    print("        Strong evidence for RH connection.")
elif n_failed <= 2:
    print("RESULT: Near-perfect discrimination. A few borderline cases warrant investigation.")
else:
    print("RESULT: Significant failures — operator may not be sufficiently selective.")
print("\nDone.")
