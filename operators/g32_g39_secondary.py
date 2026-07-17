#!/usr/bin/env python3
"""
g32_g39_secondary.py — Secondary sweeps for unconfirmed zeros.

γ₃₂ = 105.446623: p=2053 (err=0.000036) → late converger (N64=-4.9, N256=-10.3, delta=-5.36)
  Goal: test p=2053 at N=320,384,448 + search top-10 candidates with revised criterion

γ₃₉ = 121.370125: all top-4 diverge — need wider search
  Goal: search primes [50, 30000], top-8 candidates

Revised stability criterion: stable = (N256 < -10) AND (delta < +2.5)
  [one-sided — only penalize DIVERGENCE, not late convergence]
"""
import math, numpy as np, scipy.linalg
from scipy.special import digamma

TWO_PI = 2 * math.pi
LOGPI = math.log(math.pi)
G_PER_N = 900

GAMMA32 = 105.446623052
GAMMA39 = 121.370125002


def make_sieve(MAX_LSQ):
    sv = np.ones(MAX_LSQ + 1, dtype=bool)
    sv[0] = sv[1] = False
    for i in range(2, int(MAX_LSQ ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    return sv, np.where(sv)[0]


sv30000, primes30000 = make_sieve(30000)


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
    p_in = primes30000[primes30000 <= lsq]
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


def prime_candidates(gamma, p_max=30000, top_n=8):
    cands = []
    for p in primes30000:
        if p < 50 or p > p_max:
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


def run_zero(gamma, label, candidates, nvals_base=[64, 128, 192, 256], extra_n=None):
    """Test candidates with one-sided delta criterion."""
    print(f"\n{'='*60}")
    print(f"{label} = {gamma:.6f}")
    print(f"{'='*60}")
    confirmed = []
    for err, j, p, ge in candidates:
        row = {}
        for N in nvals_base:
            e = compute_ev0(p, N)
            row[N] = math.log10(abs(e)) if abs(e) > 0 else -99
        delta = row[256] - row[64]
        # ONE-SIDED: only penalize positive delta (divergence), allow negative (convergence)
        stable = (row[256] < -10) and (delta < 2.5)
        all_stable = all(row[N] < -10 for N in nvals_base)
        vals = "  ".join(f"N{N}={row[N]:>7.4f}" for N in nvals_base)
        if all_stable:
            verdict = "STABLE***"
        elif stable:
            verdict = f"STABLE(1side,Δ={delta:+.2f})"
        else:
            verdict = f"div(Δ={delta:+.2f})"
        print(f"  p={p:6d} j={j:4d} err={err:.6f}: {vals}  → {verdict}", flush=True)
        if stable and not confirmed:
            confirmed.append((j, p, err, row[256], delta))
    return confirmed


# ============================================================
# γ₃₂: Extended N test for p=2053 + wider search
# ============================================================
print("=" * 70)
print("γ₃₂ SECONDARY SWEEP")
print("=" * 70)

print("\n--- p=2053: Extended N to verify late convergence ---")
print("Prior: N64=-4.93, N128=-6.90, N192=-6.98, N256=-10.29 (delta=-5.36)")
print("N=320 max (N=512 would need ~19 GB peak; Ollama occupies 10 GB)")
for N in [320]:
    e = compute_ev0(2053, N)
    log_e = math.log10(abs(e)) if abs(e) > 0 else -99
    status = "stable" if log_e < -10 else "ABOVE"
    print(f"  N={N}: log10|ev0|={log_e:.4f}  → {status}", flush=True)

print("\n--- Top-10 candidates, one-sided criterion ---")
cands32 = prime_candidates(GAMMA32, p_max=20000, top_n=10)
for err, j, p, ge in cands32:
    print(f"  p={p:6d}, j={j:4d}, γ_eff={ge:.6f}, err={err:.6f}")

confirmed32 = run_zero(GAMMA32, "γ₃₂", cands32)
if confirmed32:
    j, p, err, log256, delta = confirmed32[0]
    print(f"\n  CONFIRMED γ₃₂: j={j}, lsq={p}, err={err:.6f}, N256={log256:.3f}, Δ={delta:+.3f}")
else:
    print("\n  γ₃₂ NOT confirmed in top-10. Trying small-prime range...")
    cands32_small = prime_candidates(GAMMA32, p_max=5000, top_n=10)
    confirmed32 = run_zero(GAMMA32, "γ₃₂ small-prime", cands32_small)


# ============================================================
# γ₃₉: Wide search up to lsq=30000, top-8 candidates
# ============================================================
print("\n\n" + "=" * 70)
print("γ₃₉ SECONDARY SWEEP")
print("=" * 70)

print("\n--- Top-8 candidates in [50, 30000], one-sided criterion ---")
cands39 = prime_candidates(GAMMA39, p_max=30000, top_n=8)
for err, j, p, ge in cands39:
    print(f"  p={p:6d}, j={j:4d}, γ_eff={ge:.6f}, err={err:.6f}")

confirmed39 = run_zero(GAMMA39, "γ₃₉", cands39)
if confirmed39:
    j, p, err, log256, delta = confirmed39[0]
    print(f"\n  CONFIRMED γ₃₉: j={j}, lsq={p}, err={err:.6f}, N256={log256:.3f}, Δ={delta:+.3f}")
else:
    print("\n  γ₃₉ NOT confirmed in [50,30000] top-8.")
    print("  Checking small-prime range explicitly...")
    cands39_small = prime_candidates(GAMMA39, p_max=5000, top_n=8)
    for err, j, p, ge in cands39_small:
        print(f"  p={p:6d}, j={j:4d}, γ_eff={ge:.6f}, err={err:.6f}")
    confirmed39 = run_zero(GAMMA39, "γ₃₉ small", cands39_small)

print("\n\nDone.")
