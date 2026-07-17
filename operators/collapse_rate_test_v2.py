#!/usr/bin/env python3
"""
collapse_rate_test_v2.py

Collapse rate test — memory-safe version.
p=797 at N=512 already confirmed: -11.2655 (STABLE).
This script continues from N=640 for Cluster B and N=512 for Cluster C.

Uses G_PER_N=600 at N>=512 to keep E matrix under ~8 GB peak.
(Standard is G_PER_N=900 but 600 is sufficient for stability trend.)
"""
import math, numpy as np, scipy.linalg, sys, time
from scipy.special import digamma

TWO_PI = 2 * math.pi
LOGPI = math.log(math.pi)

def make_sieve(MAX_LSQ):
    sv = np.ones(MAX_LSQ + 1, dtype=bool)
    sv[0] = sv[1] = False
    for i in range(2, int(MAX_LSQ**0.5) + 1):
        if sv[i]:
            sv[i*i::i] = False
    return sv, np.where(sv)[0]

sv, primes = make_sieve(30000)

def compute_ev0(lsq, N, g_per_n=600):
    t0 = time.time()
    G = g_per_n * N
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
    del E  # free ASAP after arch
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

print("Cluster Collapse Rate Test v2 (memory-safe, G_PER_N=600)", flush=True)
print(f"γ₁ = 14.1347 | Previously confirmed: p=797 N=512 = -11.2655 (STABLE)", flush=True)
print("="*70, flush=True)

# Cluster B: p=797 — skip N=512 (already done), test N=640, 768
print("\nCluster B: p=797 (j=15, nat_freq=14.108, ε=0.027 from γ₁)", flush=True)
print("  [Known: N=512 → -11.2655 STABLE]", flush=True)
print(f"  {'N':>5}  {'G_PER_N':>7}  {'log10|ev0|':>12}  {'stable':>6}  {'time':>6}", flush=True)
prev_log = -11.2655  # from prior run
for N, gpn in [(640, 600), (768, 600)]:
    ev0, elapsed = compute_ev0(797, N, g_per_n=gpn)
    log_ev0 = math.log10(abs(ev0)) if abs(ev0) > 0 else -99.0
    stable = "YES" if log_ev0 < -10 else "NO"
    delta = log_ev0 - prev_log
    print(f"  N={N:>4}  gpn={gpn}  {log_ev0:>12.4f}  {stable:>6}  {elapsed:.0f}s (Δ={delta:+.2f})", flush=True)
    prev_log = log_ev0
    if log_ev0 > -8:
        print(f"  *** DIVERGED at N={N} — cluster collapsed ***", flush=True)
        break

# Cluster C: p=19417 (ε=0.134) and p=19429 (ε=0.135) — test N=512, 640
for p, j, nat_freq, eps, label in [
    (19417, 22, 14.000, 0.134, "C-1"),
    (19429, 22, 13.999, 0.135, "C-2"),
]:
    print(f"\nCluster {label}: p={p} (j={j}, nat_freq={nat_freq:.3f}, ε={eps:.3f})", flush=True)
    print(f"  [Known: N=448 → stable from high_n_fake_test]", flush=True)
    print(f"  {'N':>5}  {'G_PER_N':>7}  {'log10|ev0|':>12}  {'stable':>6}  {'time':>6}", flush=True)
    prev_log = None
    for N, gpn in [(512, 600), (640, 600)]:
        ev0, elapsed = compute_ev0(p, N, g_per_n=gpn)
        log_ev0 = math.log10(abs(ev0)) if abs(ev0) > 0 else -99.0
        stable = "YES" if log_ev0 < -10 else "NO"
        trend = ""
        if prev_log is not None:
            delta = log_ev0 - prev_log
            trend = f" (Δ={delta:+.2f})"
        print(f"  N={N:>4}  gpn={gpn}  {log_ev0:>12.4f}  {stable:>6}  {elapsed:.0f}s{trend}", flush=True)
        prev_log = log_ev0
        if log_ev0 > -8:
            print(f"  *** DIVERGED at N={N} — cluster collapsing ***", flush=True)
            break

print("\nDone.", flush=True)
