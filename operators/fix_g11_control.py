#!/usr/bin/env python3
"""
fix_g11_control.py — Two things:

1. Fix γ₁₁ bug: use direct prime enumeration (for each prime p, compute best j)
   instead of j→lsq (roundoff can skip primes when exp(2πj/γ) rounds to composite).

2. Control experiment: test "null" primes — primes NOT near any resonant scale for
   any zero 1-20. If they're also stable, small-lsq stability is spurious. If they
   diverge, the stability IS zero-specific.

Null primes to test: 53, 71, 79, 89 — chosen so γ_eff for best j is far (>0.15)
from ALL known zeros 1-20.
"""
import math, numpy as np, scipy.linalg
from scipy.special import digamma

TWO_PI = 2 * math.pi
LOGPI = math.log(math.pi)
G_PER_N = 900

ZEROS_1_20 = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
    52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
    67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
]


def make_sieve(MAX_LSQ):
    sv = np.ones(MAX_LSQ + 1, dtype=bool)
    sv[0] = sv[1] = False
    for i in range(2, int(MAX_LSQ ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    return sv, np.where(sv)[0]


sv5000, primes5000 = make_sieve(5000)


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
    p_in = primes5000[primes5000 <= lsq]
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


def best_zero_match(lsq):
    """For a given prime lsq, find best zero match using γ_eff."""
    log_lsq = math.log(lsq)
    best_err = 999.0
    best = None
    for gam in ZEROS_1_20:
        j = round(gam * log_lsq / TWO_PI)
        ge = TWO_PI * j / log_lsq
        err = abs(ge - gam)
        if err < best_err:
            best_err = err
            best = (gam, j, ge, err)
    return best


def prime_candidates_direct(gamma, p_min=50, p_max=2500, top_n=6):
    """Enumerate primes directly, compute best j for each, rank by γ_eff error."""
    cands = []
    for p in primes5000:
        if p < p_min:
            continue
        if p > p_max:
            break
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


def test_lsq(lsq, label=""):
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
    print(f"  {label}: {vals}  → {verdict}", flush=True)
    return stable, all_stable, row


# ======================================================================
# PART 1: FIX γ₁₁ — direct prime enumeration
# ======================================================================
GAMMA11 = 52.970321477714000
print("=" * 70)
print(f"PART 1: γ₁₁ = {GAMMA11} — Direct prime enumeration")
print("(Prior j→lsq approach had roundoff bug; this enumerates primes directly)")
print("=" * 70)

cands11 = prime_candidates_direct(GAMMA11, p_min=50, p_max=2500, top_n=8)
print(f"\nTop prime candidates (direct enumeration, sorted by γ_eff error):")
for err, j, p, ge in cands11:
    print(f"  p={p:5d}, j={j:3d}, γ_eff={ge:.6f}, err={err:.6f}")

print()
g11_confirmed = None
for err, j, p, ge in cands11:
    stable, all_s, row = test_lsq(p, f"p={p:5d} j={j} err={err:.5f}")
    if stable and g11_confirmed is None:
        g11_confirmed = (p, j, ge, err, row[256], row[256] - row[64])

if g11_confirmed:
    p, j, ge, err, n256, delta = g11_confirmed
    print(f"\nγ₁₁ CONFIRMED: lsq={p}, j={j}, γ_eff={ge:.6f}, err={err:.6f}, N256={n256:.3f}")
else:
    print(f"\nγ₁₁ still not confirmed in [50, 2500] — need wider range")
    # Try p_max=5000
    print(f"\nExtending to p_max=5000...")
    cands11_wide = prime_candidates_direct(GAMMA11, p_min=50, p_max=5000, top_n=6)
    print(f"Additional candidates in [2500, 5000]:")
    for err, j, p, ge in cands11_wide:
        if p > 2500:
            print(f"  p={p:5d}, j={j:3d}, γ_eff={ge:.6f}, err={err:.6f}")
    for err, j, p, ge in cands11_wide:
        if p > 2500:
            stable, all_s, row = test_lsq(p, f"p={p:5d} j={j} err={err:.5f}")
            if stable and g11_confirmed is None:
                g11_confirmed = (p, j, ge, err, row[256], row[256] - row[64])


# ======================================================================
# PART 2: CONTROL EXPERIMENT
# Purpose: test primes NOT near any zero. If stable → spurious.
#          If diverge → stability IS zero-specific.
#
# Find "null" primes: those where min γ_eff error over all zeros 1-20 is large.
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: CONTROL EXPERIMENT — primes far from all zeros")
print("A null prime has min γ_eff error > 0.15 across all 20 zeros.")
print("If null primes are STABLE → small-lsq stability is spurious.")
print("If null primes DIVERGE → stability is zero-specific.")
print("=" * 70)

# Find null primes: primes where best zero match has err > 0.15
null_primes = []
for p in primes5000:
    if p < 50 or p > 500:
        continue
    gam, j, ge, err = best_zero_match(p)
    if err > 0.15:
        null_primes.append((err, p, gam, j, ge))

null_primes.sort(reverse=True)  # Worst match first
print(f"\nNull primes found (err > 0.15 from all 20 zeros): {len(null_primes)}")
for err, p, gam, j, ge in null_primes[:10]:
    print(f"  p={p:4d}, best_match=γ({gam:.3f}), j={j}, γ_eff={ge:.4f}, err={err:.4f}")

print(f"\nTesting top null primes at N=64..256:")
null_results = []
for err, p, gam, j, ge in null_primes[:6]:
    stable, all_s, row = test_lsq(p, f"NULL p={p:4d} (err={err:.3f} from γ={gam:.3f})")
    null_results.append((p, stable, all_s, err, row[256]))

# ======================================================================
# SYNTHESIS
# ======================================================================
print("\n" + "=" * 70)
print("SYNTHESIS")
print("=" * 70)

if g11_confirmed:
    p, j, ge, err, n256, delta = g11_confirmed
    print(f"\nγ₁₁ CONFIRMED: lsq={p}, j={j}, γ_eff={ge:.6f}, err={err:.6f}, N256={n256:.3f}")
else:
    print(f"\nγ₁₁: NOT confirmed in [50, 5000]")

stable_nulls = sum(1 for p, s, a, e, n in null_results if s)
print(f"\nControl: {stable_nulls}/{len(null_results)} null primes are stable")
if stable_nulls == 0:
    print("  → ZERO null primes stable. Stability IS zero-specific. ")
elif stable_nulls <= 1:
    print("  → Mostly not stable. Stability is predominantly zero-specific.")
else:
    print("  → MULTIPLE null primes stable. Small-lsq stability may be partially spurious.")
    print("    Need to rely on γ_eff error as primary discriminator.")

print("\nNull prime details:")
for p, stable, all_s, err, n256 in null_results:
    status = "STABLE***" if all_s else ("STABLE(N256)" if stable else "div")
    print(f"  p={p:4d}: N256={n256:.3f} → {status}")
