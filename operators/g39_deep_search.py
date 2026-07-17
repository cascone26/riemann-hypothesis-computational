#!/usr/bin/env python3
"""
g39_deep_search.py

γ₃₉ = 121.370125 — no stable prime found in [50, 30000] top-8.
Near-miss: p=2753 (err=0.002729): N64=-5.18, N128=-7.07, N192=-9.97, N256=-9.20
           Just 0.80 above threshold at N256!

Three approaches:
  A. Test p=2753 at N=320 — does it cross -10?
  B. Search primes in [30000, 80000], top-6 by γ_eff error (extended range)
  C. Test "almost-stable" candidates at extended N

Memory constraint: N≤256 for large primes (lsq≥10000), N≤320 for small primes (lsq<5000)
  N=256, D=513: E matrix = (513, 230400) × 16 bytes = 1.89 GB (safe)
  N=320, D=641: E matrix = (641, 288000) × 16 bytes = 2.96 GB (safe with 18 GB free)
"""
import math, numpy as np, scipy.linalg
from scipy.special import digamma

TWO_PI = 2 * math.pi
LOGPI = math.log(math.pi)
G_PER_N = 900
GAMMA39 = 121.370125002


def make_sieve(MAX_LSQ):
    sv = np.ones(MAX_LSQ + 1, dtype=bool)
    sv[0] = sv[1] = False
    for i in range(2, int(MAX_LSQ ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    return sv, np.where(sv)[0]


sv80000, primes80000 = make_sieve(80000)


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
    p_in = primes80000[primes80000 <= lsq]
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


# ============================================================
# A: Test p=2753 at N=320 — near-miss candidate
# ============================================================
print("=" * 70)
print(f"γ₃₉ = {GAMMA39}")
print("=" * 70)

print("\nAPPROACH A: p=2753 (near-miss) extended N test")
print("Prior: N64=-5.18, N128=-7.07, N192=-9.97, N256=-9.20")
err, j, ge = gamma_eff_error(2753, GAMMA39)
print(f"p=2753: j={j}, γ_eff={ge:.6f}, err={err:.6f}")

e = compute_ev0(2753, 320)
log_e = math.log10(abs(e)) if abs(e) > 0 else -99
print(f"  N=320: {log_e:.4f}  → {'stable' if log_e < -10 else 'ABOVE'}", flush=True)

# ============================================================
# B: Extended range search [30000, 80000]
# ============================================================
print("\nAPPROACH B: Extended prime range [30000, 80000], top-8 by γ_eff error")

cands_extended = []
for p in primes80000:
    if p < 30000 or p > 80000:
        continue
    log_p = math.log(p)
    j = round(GAMMA39 * log_p / TWO_PI)
    if j < 1:
        continue
    ge = TWO_PI * j / log_p
    err = abs(ge - GAMMA39)
    cands_extended.append((err, j, p, ge))

cands_extended.sort()
top8 = cands_extended[:8]

print(f"  Top candidates in [30000, 80000]:")
for err, j, p, ge in top8:
    print(f"    p={p:6d}, j={j:4d}, γ_eff={ge:.6f}, err={err:.6f}")

NVALS = [64, 128, 192, 256]
confirmed = []
for err, j, p, ge in top8:
    row = {}
    # Large primes at N=256 need ~1.9 GB for E — safe
    for N in NVALS:
        e = compute_ev0(p, N)
        row[N] = math.log10(abs(e)) if abs(e) > 0 else -99
    delta = row[256] - row[64]
    stable = (row[256] < -10) and (delta < 2.5)
    all_stable = all(row[N] < -10 for N in NVALS)
    vals = "  ".join(f"N{N}={row[N]:>7.4f}" for N in NVALS)
    verdict = "STABLE***" if all_stable else ("STABLE(1side)" if stable else f"div(Δ={delta:+.2f})")
    print(f"  p={p:6d} j={j:4d} err={err:.6f}: {vals}  → {verdict}", flush=True)
    if stable and not confirmed:
        confirmed.append((j, p, err, ge, row[256], delta))

# ============================================================
# C: Revisit [50, 30000] with wider error threshold (err < 0.005)
# ============================================================
print("\nAPPROACH C: Wider error threshold in [50, 30000], top-12 including err<0.005")

cands_wide = []
for p in primes80000:
    if p < 50 or p > 30000:
        continue
    log_p = math.log(p)
    j = round(GAMMA39 * log_p / TWO_PI)
    if j < 1:
        continue
    ge = TWO_PI * j / log_p
    err = abs(ge - GAMMA39)
    if err < 0.005:
        cands_wide.append((err, j, p, ge))

cands_wide.sort()
new_cands = [(err, j, p, ge) for err, j, p, ge in cands_wide if err >= 0.001][:8]
print(f"  New candidates with 0.001 ≤ err < 0.005 in [50, 30000]:")
for err, j, p, ge in new_cands:
    print(f"    p={p:6d}, j={j:4d}, γ_eff={ge:.6f}, err={err:.6f}")

for err, j, p, ge in new_cands:
    row = {}
    n_max = 320 if p < 5000 else 256  # memory limit
    nvals = [64, 128, 192, n_max]
    for N in nvals:
        e = compute_ev0(p, N)
        row[N] = math.log10(abs(e)) if abs(e) > 0 else -99
    delta = row[n_max] - row[64]
    stable = (row[n_max] < -10) and (delta < 2.5)
    all_stable = all(row[N] < -10 for N in nvals)
    vals = "  ".join(f"N{N}={row[N]:>7.4f}" for N in nvals)
    verdict = "STABLE***" if all_stable else ("STABLE(1side)" if stable else f"div(Δ={delta:+.2f})")
    print(f"  p={p:6d} j={j:4d} err={err:.6f}: {vals}  → {verdict}", flush=True)
    if stable and not confirmed:
        confirmed.append((j, p, err, ge, row[n_max], delta))

print("\n")
if confirmed:
    j, p, err, ge, logN, delta = confirmed[0]
    print(f"CONFIRMED γ₃₉: j={j}, lsq={p}, err={err:.6f}, ev₀(N256/320)={logN:.3f}, Δ={delta:+.3f}")
else:
    print("γ₃₉ NOT confirmed in extended search.")
    print("Near-miss: p=2753 (err=0.002729) reached N192=-9.97 (0.03 below threshold)")
    print("May need even broader range or this zero has unusual structure.")

print("\nDone.")
