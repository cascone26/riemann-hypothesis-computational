#!/usr/bin/env python3
"""
zeros_31_40_scan.py — Detect Riemann zeros γ₃₁..γ₄₀.

Values from Odlyzko tables (9+ decimal places).
Method: direct prime enumeration, top-4 by γ_eff error, N=64..256.
"""
import math, numpy as np, scipy.linalg
from scipy.special import digamma

TWO_PI = 2 * math.pi
LOGPI = math.log(math.pi)
G_PER_N = 900

ZEROS_31_40 = [
    (31, 103.725538040),
    (32, 105.446623052),
    (33, 107.168611184),
    (34, 111.029535543),
    (35, 111.874659178),
    (36, 114.320220915),
    (37, 116.226680321),
    (38, 118.790782866),
    (39, 121.370125002),
    (40, 122.946829294),
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


def prime_candidates(gamma, p_min=50, p_max=20000, top_n=4):
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


NVALS = [64, 128, 192, 256]
confirmed = []

print("=" * 70)
print("CvS Operator Scan: Riemann Zeros γ₃₁..γ₄₀")
print("Method: direct prime enumeration, top 4 by γ_eff error")
print("Stability: log(N=256) < -10 AND |Δ(64→256)| < 2.5")
print("=" * 70)

for k, gamma in ZEROS_31_40:
    print(f"\n{'='*60}")
    print(f"γ{k} = {gamma:.6f}")
    print(f"{'='*60}")

    cands = prime_candidates(gamma, p_min=50, p_max=20000, top_n=4)
    print(f"  Top prime candidates:")
    for err, j, p, ge in cands:
        print(f"    p={p:6d}, j={j:4d}, γ_eff={ge:.6f}, err={err:.6f}")

    print()
    zero_confirmed = False
    for err, j, p, ge in cands:
        row = {}
        for N in NVALS:
            e = compute_ev0(p, N)
            log_e = math.log10(abs(e)) if abs(e) > 0 else -99
            row[N] = log_e
        delta = row[256] - row[64]
        all_stable = all(row[N] < -10 for N in NVALS)
        stable = (row[256] < -10) and (abs(delta) < 2.5)
        vals = "  ".join(f"N{N}={row[N]:>7.4f}" for N in NVALS)
        verdict = "STABLE***" if all_stable else ("STABLE(N256)" if stable else f"div(Δ={delta:+.2f})")
        print(f"  p={p:6d} j={j:4d} err={err:.6f}: {vals}  → {verdict}", flush=True)
        if stable and not zero_confirmed:
            zero_confirmed = True
            confirmed.append((k, gamma, j, p, err, ge, row[256], delta))

print("\n\n" + "=" * 70)
print("SYNTHESIS")
print("=" * 70)
print(f"\nZeros γ₃₁..γ₄₀: {len(confirmed)}/10 confirmed")
for k, gamma, j, p, err, ge, log256, delta in confirmed:
    print(f"  γ{k} ({gamma:.4f}): j={j}, lsq={p}, err={err:.6f}, N256={log256:.3f}, Δ={delta:+.3f}")
missing = [k for k, _ in ZEROS_31_40 if k not in [c[0] for c in confirmed]]
if missing: print(f"\nNot confirmed: {missing}")
print("\nDone.")
