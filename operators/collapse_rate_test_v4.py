#!/usr/bin/env python3
"""
collapse_rate_test_v4.py — streaming arch (never stores full E matrix)

Memory improvement over v3: combine E-build and arch-accumulate into ONE
chunk loop. E is never fully allocated. Peak = ~4 GB instead of 31 GB,
enabling N=1024+ with 22 GB available RAM.

Tests:
  TRUE resonance:  p=787  (γ₁=14.135, j=15, nat_freq=14.134, err=0.001)  → should stay stable/deeper
  Fake cluster B:  p=797  (j=15, nat_freq=14.108, ε=0.027 from γ₁)       → mystery stable prime
  Fake cluster C:  p=19417 (j=22, nat_freq=14.000, ε=0.134)              → should eventually diverge
  Fake cluster C2: p=19429 (j=22, nat_freq=13.999, ε=0.135)              → should eventually diverge

Prior results (v3):
  p=787:  (to test)
  p=797:  N=512 -11.2655  N=640 -11.5628  N=768 -10.5629
  p=19417: N=448 stable  N=512 -10.3608  N=640 -10.5156
  p=19429: N=448 stable  N=512 -10.7201  N=640 -10.3211
"""
import math, numpy as np, scipy.linalg, sys, time
from scipy.special import digamma

TWO_PI = 2 * math.pi
LOGPI = math.log(math.pi)
G_PER_N = 900
CHUNK_COLS = 24000  # smaller chunks to reduce peak (3 temps × chunk at once)

def make_sieve(MAX_LSQ):
    sv = np.ones(MAX_LSQ + 1, dtype=bool)
    sv[0] = sv[1] = False
    for i in range(2, int(MAX_LSQ**0.5)+1):
        if sv[i]:
            sv[i*i::i] = False
    return np.where(sv)[0]

primes = make_sieve(30000)

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
    wOm = w * Om

    c = TWO_PI * ns / L  # (D,)

    # Streaming arch: build E chunk-by-chunk, accumulate arch WITHOUT storing full E
    arch = np.zeros((D, D), dtype=float)
    for g0 in range(0, G, CHUNK_COLS):
        g1 = min(g0 + CHUNK_COLS, G)
        d = c[:, None] - grid[None, g0:g1]          # D × chunk, real
        sm = np.abs(d) < 1e-12
        d_safe = np.where(sm, 1.0, d)               # avoid /0
        Ec = np.where(
            sm,
            np.sqrt(L) + 0j,
            (1/np.sqrt(L)) * np.exp(0.5j * d * L) * 2 * np.sin(0.5 * d * L) / d_safe
        )  # D × chunk, complex
        Ew = Ec * wOm[g0:g1]
        arch += (Ew @ Ec.conj().T).real / TWO_PI
        del d, sm, d_safe, Ec, Ew

    sh = math.sinh(L / 4) ** 2
    ni = ns[:, None]; nj = ns[None, :]
    B = 32 * L * sh * (L**2 - 16*math.pi**2 * ni*nj) / (
        (L**2 + 16*math.pi**2 * ni**2) * (L**2 + 16*math.pi**2 * nj**2))

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
        S  += np.sin(ph) @ wts_p[k0:k1]
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

def run_cluster(name, p, j, nat_freq, eps, desc, N_list, prior_N, prior_log):
    print(f"\n{name}: p={p} (j={j}, nat_freq={nat_freq:.3f}, {desc})", flush=True)
    if prior_N:
        print(f"  Prior: N={prior_N} → {prior_log:.4f} STABLE", flush=True)
    print(f"  {'N':>6}  {'log10|ev0|':>12}  {'stable':>6}  {'time':>6}  {'delta':>8}", flush=True)
    prev_log = prior_log
    for N in N_list:
        ev0, elapsed = compute_ev0(p, N)
        log_ev0 = math.log10(abs(ev0)) if abs(ev0) > 0 else -99.0
        stable = "YES" if log_ev0 < -10 else "NO"
        delta_str = f"Δ={log_ev0 - prev_log:+.2f}" if prev_log is not None else ""
        print(f"  N={N:>5}: {log_ev0:>12.4f}  {stable:>6}  {elapsed:.0f}s  {delta_str}", flush=True)
        sys.stdout.flush()
        prev_log = log_ev0
        if log_ev0 > -8:
            print(f"  *** DIVERGED at N={N} — clearly non-resonant ***", flush=True)
            break

print("Collapse Rate Test v4 (streaming arch, CHUNK_COLS=24000)", flush=True)
print(f"G_PER_N={G_PER_N} | Streaming: arch computed in same chunk loop as E, E never fully stored", flush=True)
print(f"Peak memory estimate: ~4 GB (vs ~31 GB for v3 at N=1024)", flush=True)
print("="*70, flush=True)

# True resonance: p=787, the confirmed stable prime for γ₁=14.135
run_cluster(
    "TRUE RESONANCE", 787, 15, 14.134, 0.001,
    "ε=0.001 from γ₁=14.135 (CONFIRMED STABLE PRIME)",
    [768, 896, 1024, 1280],
    512, -11.5  # prior approx from earlier sessions
)

# Cluster B: p=797 (mystery prime, stable but not a true resonance)
run_cluster(
    "Cluster B (mystery)", 797, 15, 14.108, 0.027,
    "ε=0.027 from γ₁=14.135",
    [896, 1024, 1280],
    768, -10.5629
)

# Cluster C-1: p=19417 (fake stable near γ_fake=14.0)
run_cluster(
    "Cluster C-1 (fake)", 19417, 22, 14.000, 0.134,
    "ε=0.134 from γ₁=14.135, nat_freq of fake γ=14.0",
    [768, 896, 1024],
    640, -10.5156
)

# Cluster C-2: p=19429 (fake stable near γ_fake=14.0)
run_cluster(
    "Cluster C-2 (fake)", 19429, 22, 13.999, 0.135,
    "ε=0.135 from γ₁=14.135",
    [768, 896, 1024],
    640, -10.3211
)

print("\nDone.", flush=True)
