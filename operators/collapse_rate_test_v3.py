#!/usr/bin/env python3
"""
collapse_rate_test_v3.py

Memory-efficient version using chunked E construction + chunked arch.
Peak memory: ~E + chunk (~E/8) instead of 4×E.

Previously confirmed:
  p=797 N=512 = -11.2655 (STABLE)  [G_PER_N=900]

Tests:
  Cluster B: p=797 at N=640, 768 (G_PER_N=900)
  Cluster C: p=19417, p=19429 at N=512, 640 (G_PER_N=900)
"""
import math, numpy as np, scipy.linalg, sys, time
from scipy.special import digamma

TWO_PI = 2 * math.pi
LOGPI = math.log(math.pi)
G_PER_N = 900
CHUNK_COLS = 48000  # columns per chunk (reduces peak memory from 4×E to ~E + E/8)

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
    wOm = w * Om  # shape (G,) — quadrature weights × Omega

    # Build E in-place chunk by chunk to keep peak at ~E + chunk
    c = TWO_PI * ns / L  # shape (D,)
    E = np.empty((D, G), dtype=complex)
    for g0 in range(0, G, CHUNK_COLS):
        g1 = min(g0 + CHUNK_COLS, G)
        d_chunk = c[:, None] - grid[None, g0:g1]  # D × chunk
        sm = np.abs(d_chunk) < 1e-12
        rhs = np.where(sm, np.sqrt(L) + 0j,
                       (1/np.sqrt(L)) * np.exp(0.5j * d_chunk * L) *
                       2 * np.sin(0.5 * d_chunk * L) / np.where(sm, 1.0, d_chunk))
        E[:, g0:g1] = rhs
    # E is now fully built. Peak so far: E (full) + d_chunk (~E/8) + temp (~E/8)

    # Compute arch in chunks: arch += (E[:,chunk] * wOm[chunk]) @ E[:,chunk].H
    arch = np.zeros((D, D), dtype=float)
    for g0 in range(0, G, CHUNK_COLS):
        g1 = min(g0 + CHUNK_COLS, G)
        Ec = E[:, g0:g1]                    # D × chunk (view, no copy)
        Ew = Ec * wOm[g0:g1]               # D × chunk (new array ~E/8)
        arch += (Ew @ Ec.conj().T).real / TWO_PI  # D × D (tiny)
    del E, Ew, Ec

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

print("Cluster Collapse Rate Test v3 (chunked, G_PER_N=900)", flush=True)
print(f"γ₁ = 14.1347 | CHUNK_COLS={CHUNK_COLS}", flush=True)
print(f"Known: p=797 N=512 = -11.2655 (STABLE) [G_PER_N=900]", flush=True)
print("="*70, flush=True)

# Cluster B: p=797 — start at N=640 (N=512 already done)
print("\nCluster B: p=797 (j=15, nat_freq=14.108, ε=0.027 from γ₁=14.135)", flush=True)
print("  [N=512 → -11.2655 STABLE (confirmed)]", flush=True)
print(f"  {'N':>5}  {'log10|ev0|':>12}  {'stable':>6}  {'time':>6}  {'delta':>8}", flush=True)
prev_log = -11.2655
for N in [640, 768]:
    ev0, elapsed = compute_ev0(797, N)
    log_ev0 = math.log10(abs(ev0)) if abs(ev0) > 0 else -99.0
    stable = "YES" if log_ev0 < -10 else "NO"
    delta = log_ev0 - prev_log
    print(f"  N={N:>4}: {log_ev0:>12.4f}  {stable:>6}  {elapsed:.0f}s  Δ={delta:+.2f}", flush=True)
    prev_log = log_ev0
    if log_ev0 > -8:
        print(f"  *** DIVERGED at N={N} — cluster collapsed ***", flush=True)
        break

# Cluster C: p=19417 and p=19429 — N=512, 640
for p, j, nat_freq, eps, label in [
    (19417, 22, 14.000, 0.134, "C-1"),
    (19429, 22, 13.999, 0.135, "C-2"),
]:
    print(f"\nCluster {label}: p={p} (j={j}, nat_freq={nat_freq:.3f}, ε={eps:.3f} from γ_fake=14.0)", flush=True)
    print(f"  [N=448 → stable from high_n_fake_test]", flush=True)
    print(f"  {'N':>5}  {'log10|ev0|':>12}  {'stable':>6}  {'time':>6}  {'delta':>8}", flush=True)
    prev_log = None
    for N in [512, 640]:
        ev0, elapsed = compute_ev0(p, N)
        log_ev0 = math.log10(abs(ev0)) if abs(ev0) > 0 else -99.0
        stable = "YES" if log_ev0 < -10 else "NO"
        trend = ""
        if prev_log is not None:
            delta = log_ev0 - prev_log
            trend = f"  Δ={delta:+.2f}"
        print(f"  N={N:>4}: {log_ev0:>12.4f}  {stable:>6}  {elapsed:.0f}s{trend}", flush=True)
        prev_log = log_ev0
        if log_ev0 > -8:
            print(f"  *** DIVERGED at N={N} — cluster collapsing ***", flush=True)
            break

print("\nDone.", flush=True)
