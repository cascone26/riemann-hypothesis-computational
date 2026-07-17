#!/usr/bin/env python3
"""
small_lsq_scan_11_20.py — Find true resonances for the 6 failing zeros from zeros_11_20_scan.

INSIGHT: For zeros 1-10, true resonances were all at small lsq (31-491).
Large-lsq candidates are noisy (many primes = noisy P-matrix).

For zero γ_k, lsq ∈ [50, 2000] corresponds to j ∈ [γ_k×log(50)/(2π), γ_k×log(2000)/(2π)].
We enumerate ALL j in that range, find prime lsq candidates, rank by γ_eff error.

Failing zeros: γ₁₁, γ₁₂, γ₁₅, γ₁₆, γ₁₉, γ₂₀
"""
import math, numpy as np, scipy.linalg
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
    return sv, np.where(sv)[0]


sv3000, primes3000 = make_sieve(3000)


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
    p_in = primes3000[primes3000 <= lsq]
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


def find_prime_candidates(gamma, lsq_min=50, lsq_max=2500, top_n=6):
    """Find top prime lsq candidates by γ_eff error."""
    j_min = int(gamma * math.log(lsq_min) / TWO_PI) - 1
    j_max = int(gamma * math.log(lsq_max) / TWO_PI) + 2
    candidates = []
    seen_lsq = set()
    for j in range(max(1, j_min), j_max + 1):
        lsq = round(math.exp(TWO_PI * j / gamma))
        if lsq < lsq_min or lsq > lsq_max:
            continue
        if lsq in seen_lsq or lsq > 3000:
            continue
        seen_lsq.add(lsq)
        if not sv3000[lsq]:
            continue
        g_eff = TWO_PI * j / math.log(lsq)
        err = abs(g_eff - gamma)
        candidates.append((err, j, lsq, g_eff))
    candidates.sort()
    return candidates[:top_n]


FAILING = [
    (11, 52.970321477714000),
    (12, 56.446247697086000),
    (15, 65.112544048081000),
    (16, 67.079810529494000),
    (19, 75.704690699083000),
    (20, 77.144840068875000),
]

NVALS = [64, 128, 192, 256]
confirmed = []

print("=" * 70)
print("Small-lsq scan for 6 failing zeros from zeros_11_20_scan")
print("lsq range: [50, 2500] | top 6 prime candidates by γ_eff error")
print("=" * 70)

for k, gamma in FAILING:
    print(f"\n{'='*60}")
    print(f"γ₁{k-10} = {gamma:.6f}")
    print(f"{'='*60}")

    cands = find_prime_candidates(gamma, lsq_min=50, lsq_max=2500, top_n=6)
    if not cands:
        print("  No prime candidates found in range [50, 2500]!")
        continue

    print(f"  Top candidates (prime, sorted by γ_eff error):")
    for err, j, lsq, ge in cands:
        print(f"    j={j}, lsq={lsq}, γ_eff={ge:.6f}, err={err:.6f}")

    print()
    zero_confirmed = False
    for err, j, lsq, ge in cands:
        row = {}
        for N in NVALS:
            e = compute_ev0(lsq, N)
            log_e = math.log10(abs(e)) if abs(e) > 0 else -99
            row[N] = log_e
        delta = row[256] - row[64]
        all_stable = all(row[N] < -10 for N in NVALS)
        stable = (row[256] < -10) and (abs(delta) < 2.5)
        vals = "  ".join(f"N{N}={row[N]:>7.4f}" for N in NVALS)
        verdict = "STABLE***" if all_stable else ("STABLE(N256)" if stable else f"div(Δ={delta:+.2f})")
        print(f"  j={j:3d} lsq={lsq:5d} err={err:.5f}: {vals}  → {verdict}", flush=True)
        if stable and not zero_confirmed:
            zero_confirmed = True
            confirmed.append((k, gamma, j, lsq, err, ge, row[256], delta))

print("\n\n" + "=" * 70)
print("SYNTHESIS")
print("=" * 70)
print(f"\nConfirmed from small-lsq scan ({len(confirmed)}/6):")
for k, gamma, j, lsq, err, ge, log256, delta in confirmed:
    print(f"  γ₁{k-10} ({gamma:.4f}): j={j}, lsq={lsq}, γ_eff={ge:.6f}, err={err:.6f}, N256={log256:.3f}, Δ={delta:+.3f}")

still_missing = [k for k, _ in FAILING if k not in [c[0] for c in confirmed]]
if still_missing:
    print(f"\nStill not confirmed: γ indices {still_missing}")
    print("  → Need wider lsq range or larger j search")

print("\nDone.")
