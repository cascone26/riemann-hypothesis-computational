#!/usr/bin/env python3
"""
collapse_rate_test.py

Tests the Cluster Collapse Rate Law:
  HYPOTHESIS: N_collapse(cluster) ∝ 1/ε where ε = |nat_freq - nearest_zero|

Test cases:
  Cluster A: p≈787, ε≈0.001 from γ₁=14.135 → collapsed by N=256 (known)
  Cluster B: p=797, ε≈0.027 from γ₁=14.135 → stable through N=512 (mystery)
  Cluster C: p≈19417, ε≈0.134 from γ₁=14.135 → two stable at N=448

Tests C cluster at N=512, 640 to determine collapse N.
Also tests B (p=797) at N=512, 640, 768 to see if it diverges.
"""
import math, numpy as np, scipy.linalg, sys, time
from scipy.special import digamma

TWO_PI = 2 * math.pi
LOGPI = math.log(math.pi)
G_PER_N = 900

def make_sieve(MAX_LSQ):
    sv = np.ones(MAX_LSQ + 1, dtype=bool)
    sv[0] = sv[1] = False
    for i in range(2, int(MAX_LSQ**0.5) + 1):
        if sv[i]:
            sv[i*i::i] = False
    return sv, np.where(sv)[0]

sv, primes = make_sieve(30000)

def compute_ev0(lsq, N):
    t0 = time.time()
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
    sh = math.sinh(L / 4) ** 2
    ni = ns[:, None]; nj = ns[None, :]
    B = 32 * L * sh * (L**2 - 16 * math.pi**2 * ni * nj) / (
        (L**2 + 16 * math.pi**2 * ni**2) * (L**2 + 16 * math.pi**2 * nj**2))
    p_in = primes[primes <= lsq]
    ll, ww = [], []
    for p in p_in:
        lp = math.log(p); pk = p
        while pk <= lsq:
            ll.append(math.log(pk)); ww.append(lp / math.sqrt(pk))
            pk *= p
    lkpp = np.array(ll); wts_p = np.array(ww); K = len(lkpp)
    S = np.zeros(D); Sd = np.zeros(D)
    for k0 in range(0, K, 10000):
        k1 = min(k0 + 10000, K)
        ph = TWO_PI * np.outer(ns, lkpp[k0:k1]) / L
        S += np.sin(ph) @ wts_p[k0:k1]
        Sd += np.cos(ph) @ (2 * wts_p[k0:k1] * (1 - lkpp[k0:k1] / L))
    dn = ns[:, None] - ns[None, :]
    Si = S[:, None]; Sj = S[None, :]
    with np.errstate(divide='ignore', invalid='ignore'):
        P = np.where(dn != 0, (Sj - Si) / (math.pi * dn), 0.0)
    np.fill_diagonal(P, Sd)
    tau = (arch + B - P)
    tau = (tau + tau.T) / 2
    ev0 = float(scipy.linalg.eigvalsh(tau, subset_by_index=[0, 0])[0])
    return ev0, time.time() - t0

print("Cluster Collapse Rate Test")
print(f"γ₁ = 14.1347 | G_PER_N={G_PER_N}")
print("="*70)

# Test cases: (prime, nat_freq_j15_or_j22, label, epsilon)
CLUSTERS = [
    # Cluster B: mystery prime p=797 (ε=0.027 from γ₁)
    # At j=15: 2π×15/log(797) = 14.108
    [(797, 15, 14.108, 0.027, "B (mystery)")],
    # Cluster C: fake-stable primes for γ_fake=14.0 (ε≈0.134 from γ₁)
    # At j=22: nat_freq ≈ 14.000
    [(19417, 22, 14.000, 0.134, "C-1 (fake-stable)"),
     (19429, 22, 13.999, 0.135, "C-2 (fake-stable)")],
]

N_VALUES_B = [512, 640, 768]    # for cluster B (smaller prime, faster)
N_VALUES_C = [512, 640]          # for cluster C (large prime, slower)

for cluster_list in CLUSTERS:
    for p, j, nat_freq, eps, label in cluster_list:
        print(f"\np={p} (j={j}, nat_freq={nat_freq:.3f}, ε={eps:.3f}) — Cluster {label}")
        print(f"  {'N':>5}  {'log10|ev0|':>12}  {'stable':>6}  {'time':>6}")
        prev_log = None
        nvals = N_VALUES_B if p < 1000 else N_VALUES_C
        for N in nvals:
            ev0, elapsed = compute_ev0(p, N)
            log_ev0 = math.log10(abs(ev0)) if abs(ev0) > 0 else -99.0
            stable = "YES" if log_ev0 < -10 else "NO"
            trend = ""
            if prev_log is not None:
                delta = log_ev0 - prev_log
                trend = f" (Δ={delta:+.2f})"
            print(f"  N={N:>4}: {log_ev0:>12.4f}  {stable:>6}  {elapsed:.0f}s{trend}", flush=True)
            prev_log = log_ev0
            if log_ev0 > -10:
                print(f"  *** DIVERGED at N={N} ***", flush=True)
                break

print("\nDone.")
