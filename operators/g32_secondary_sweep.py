#!/usr/bin/env python3
"""
g32_secondary_sweep.py

γ₃₂ = 105.446623 — top-4 candidates from main scan all failed.
Best: p=2053, err=0.000036 → N256=-10.29, but started at N64=-4.9 (late converger).
      delta=-5.36 (converging) → wrongly classified as div by abs(delta)<2.5 criterion.

Two goals:
  1. Test p=2053 at N=320,384,448 — is the convergence real or will it revert?
  2. Expand to top-12 candidates — find any other stable prime for γ₃₂
  3. Restate stability criterion: only penalize POSITIVE delta (diverging), not negative (converging)
     Revised: stable = (N256 < -10) AND (delta < 2.5)  [one-sided]
"""
import math, numpy as np, scipy.linalg
from scipy.special import digamma

TWO_PI = 2 * math.pi
LOGPI = math.log(math.pi)
G_PER_N = 900
GAMMA32 = 105.446623052


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


def prime_candidates(gamma, top_n=12):
    cands = []
    for p in primes20000:
        if p < 50:
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


print("=" * 70)
print(f"γ₃₂ = {GAMMA32}")
print("Secondary sweep: top-12 candidates + extended N for p=2053")
print("=" * 70)

print("\n--- PART 1: Extended N sweep for p=2053 (late converger) ---")
print("Prior: N64=-4.93, N128=-6.90, N192=-6.98, N256=-10.29")
print("Testing N=320,384,448 — does convergence hold?")
for N in [320, 384, 448]:
    e = compute_ev0(2053, N)
    log_e = math.log10(abs(e)) if abs(e) > 0 else -99
    status = "stable" if log_e < -10 else "ABOVE"
    print(f"  N={N}: {log_e:.4f}  → {status}", flush=True)

print("\n--- PART 2: Top-12 candidates, revised stability criterion ---")
print("Revised: stable = (N256 < -10) AND (delta < 2.5) [one-sided — only penalize rising]")

cands = prime_candidates(GAMMA32, top_n=12)
NVALS = [64, 128, 192, 256]
confirmed = []

for err, j, p, ge in cands:
    row = {}
    for N in NVALS:
        e = compute_ev0(p, N)
        log_e = math.log10(abs(e)) if abs(e) > 0 else -99
        row[N] = log_e
    delta = row[256] - row[64]
    # Revised criterion: only penalize POSITIVE delta (divergence)
    stable_revised = (row[256] < -10) and (delta < 2.5)
    all_stable = all(row[N] < -10 for N in NVALS)
    vals = "  ".join(f"N{N}={row[N]:>7.4f}" for N in NVALS)
    if all_stable:
        verdict = "STABLE***"
    elif stable_revised:
        verdict = f"STABLE(revised,Δ={delta:+.2f})"
    else:
        verdict = f"div(Δ={delta:+.2f})"
    print(f"  p={p:6d} j={j:4d} err={err:.6f}: {vals}  → {verdict}", flush=True)
    if stable_revised and not confirmed:
        confirmed.append((j, p, err, ge, row[256], delta))

print()
if confirmed:
    j, p, err, ge, log256, delta = confirmed[0]
    print(f"CONFIRMED: j={j}, lsq={p}, err={err:.6f}, N256={log256:.3f}, Δ={delta:+.3f}")
else:
    print("Still no stable candidate in top-12. γ₃₂ needs wider range or higher N.")
print("\nDone.")
