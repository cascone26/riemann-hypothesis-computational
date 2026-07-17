#!/usr/bin/env python3
"""
zeros_11_20_scan.py — Test prime-lsq resonances for γ₁₁..γ₂₀.

Best prime candidates found by minimizing |γ_eff - γ_k| within lsq < 20000:
  γ₁₁ = 52.970321: j=78,  lsq=10427, err=0.000120
  γ₁₂ = 56.446248: j=76,  lsq=4721,  err=0.000062
  γ₁₃ = 59.347044: j=83,  lsq=6551,  err=0.000026
  γ₁₄ = 60.831779: j=86,  lsq=7207,  err=0.000347
  γ₁₅ = 65.112544: j=97,  lsq=11617, err=0.000099
  γ₁₆ = 67.079811: j=101, lsq=12841, err=0.000003
  γ₁₇ = 69.546402: j=107, lsq=15787, err=0.000022
  γ₁₈ = 72.067158: j=90,  lsq=2557,  err=0.000668
  γ₁₉ = 75.704691: j=105, lsq=6091,  err=0.000094
  γ₂₀ = 77.144840: j=98,  lsq=2927,  err=0.000329

All errors < 0.001 — much better than some confirmed zeros 1-10.

Stability criterion: log₁₀|ev₀(N=256)| < -10 AND |Δlog(64→256)| < 2.5
"""
import math, numpy as np, scipy.linalg, sys
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


# Single large sieve covering all lsq values
primes_all = make_sieve(20000)


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
    p_in = primes_all[primes_all <= lsq]
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


SCAN_TARGETS = [
    (11, 52.970321477714000, 78,  10427),
    (12, 56.446247697086000, 76,  4721),
    (13, 59.347044002602000, 83,  6551),
    (14, 60.831778524609000, 86,  7207),
    (15, 65.112544048081000, 97,  11617),
    (16, 67.079810529494000, 101, 12841),
    (17, 69.546401711174000, 107, 15787),
    (18, 72.067157674481000, 90,  2557),
    (19, 75.704690699083000, 105, 6091),
    (20, 77.144840068875000, 98,  2927),
]

NVALS = [64, 128, 192, 256]
results = {}

print("=" * 70)
print("CvS Operator Scan: Riemann Zeros γ₁₁..γ₂₀")
print("Stability: log(N=256) < -10 AND |Δ(64→256)| < 2.5")
print("=" * 70)

stable_count = 0
confirmed = []

for k, gam, j, lsq in SCAN_TARGETS:
    ge = TWO_PI * j / math.log(lsq)
    err = abs(ge - gam)
    print(f"\n{'='*60}")
    print(f"γ₁{k-10} = {gam:.6f}  (j={j}, lsq={lsq}, γ_eff={ge:.6f}, err={err:.6f})")
    print(f"{'='*60}")

    row = {}
    for N in NVALS:
        e = compute_ev0(lsq, N)
        log_e = math.log10(abs(e)) if abs(e) > 0 else -99
        row[N] = log_e
        print(f"  N={N:>4}: {log_e:>8.4f}", flush=True)

    delta = row[256] - row[64]
    stable = (row[256] < -10) and (abs(delta) < 2.5)
    # Also check all intermediate values
    all_stable = all(row[N] < -10 for N in NVALS)

    verdict = "STABLE (ALL N)" if all_stable else (
        "STABLE (N=256)" if stable else f"DIVERGING (Δ={delta:+.2f}, N256={row[256]:.2f})")
    print(f"  → {verdict}")

    results[f"g{k}_lsq{lsq}"] = {"k": k, "gamma": gam, "j": j, "lsq": lsq, "data": row,
                                    "stable": stable, "all_stable": all_stable}
    if stable:
        stable_count += 1
        confirmed.append((k, gam, j, lsq, delta, row[256]))

print("\n\n" + "=" * 70)
print("=== SYNTHESIS ===")
print("=" * 70)
print(f"\nZeros γ₁₁..γ₂₀: {stable_count}/10 show stability at N=256")
print("\nConfirmed (stable):")
for k, gam, j, lsq, delta, log256 in confirmed:
    ge = TWO_PI * j / math.log(lsq)
    err = abs(ge - gam)
    print(f"  γ₁{k-10} ({gam:.4f}): j={j}, lsq={lsq}, N256={log256:.3f}, Δ={delta:+.3f}, err={err:.6f}")

not_confirmed = [(k, gam, j, lsq) for k, gam, j, lsq in SCAN_TARGETS
                 if k not in [c[0] for c in confirmed]]
if not_confirmed:
    print("\nNot confirmed:")
    for k, gam, j, lsq in not_confirmed:
        ge = TWO_PI * j / math.log(lsq)
        r = results[f"g{k}_lsq{lsq}"]["data"]
        print(f"  γ₁{k-10} ({gam:.4f}): j={j}, lsq={lsq}, N64={r[64]:.2f}, N256={r[256]:.2f}")

import json, os
os.makedirs("results", exist_ok=True)
with open("results/zeros_11_20_scan.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nSaved → results/zeros_11_20_scan.json")
