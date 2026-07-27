#!/usr/bin/env python3
"""§356: the k=3 rung — measurements against the §355 pre-registered criteria (unchanged).

A3: slope of [S3(H) - H^2 + corrections]/H ... operationally: fit d(S3 excess)/d(log H) where
    excess(H) = [S3(H) - S3_leading(H)] / H with S3_leading = the H^2 leading term measured by
    the fit itself (fit S3(H) = a*H^2 + b*H*log H + c*H; report b). Pairwise prediction b=-1.5.
B3: same fit on S3 - S3pair; report its b (pairwise-exhausted: |b|<0.1).
C3: skewness of psi(x+H)-psi(x)-H.
S2(h) = 2 C2 prod_{p|h,p>2}(p-1)/(p-2) (h even), S3 via nu(p) on {0,d1,d2}.
"""
import numpy as np

HMAX = 4000
print(f"sieving spf to {HMAX} and building S2 ...", flush=True)
spf = np.zeros(HMAX + 1, dtype=np.int64)
for p in range(2, int(HMAX ** 0.5) + 1):
    if spf[p] == 0:
        sl = spf[p * p::p]
        sl[sl == 0] = p
        spf[p * p::p] = sl

def odd_prime_factors(m):
    fs = []
    while m % 2 == 0:
        m //= 2
    while m > 1:
        p = spf[m] if spf[m] != 0 else m
        fs.append(int(p))
        while m % p == 0:
            m //= p
    return fs

# C2 and C3 constants from primes up to 1e6
print("computing C2, C3 tail constants ...", flush=True)
sieve = np.ones(1_000_001, dtype=bool)
sieve[:2] = False
for p in range(2, 1001):
    if sieve[p]:
        sieve[p * p::p] = False
primes = np.nonzero(sieve)[0]
op = primes[primes > 2].astype(np.float64)
C2 = float(np.prod(1.0 - 1.0 / (op - 1.0) ** 2))
# S3 generic factor over p>=5 (p=3 handled per-pair: its generic factor is 0 = admissibility)
op5 = op[op >= 5]
C3 = float(np.prod((1.0 - 3.0 / op5) / (1.0 - 1.0 / op5) ** 3))
# p=2 factor for S3 with d1,d2 even: nu(2)=1: (1-1/2)/(1-1/2)^3 = 4
print(f"  C2 = {C2:.10f}   C3 (p>=5 generic product) = {C3:.10f}")

def S2(h):
    if h % 2 != 0 or h == 0:
        return 0.0
    prod = 2.0 * C2
    for p in odd_prime_factors(h):
        prod *= (p - 1.0) / (p - 2.0)
    return prod

def S3(d1, d2):
    """d1,d2 even, distinct, nonzero. p=2: nu=1 -> factor 4. p=3: actual factor directly
    (generic would be 0 -> admissibility: if {0,d1,d2} covers all residues mod 3, S3=0)."""
    r1, r2 = d1 % 3, d2 % 3
    nu3 = len({0, r1, r2})
    if nu3 == 3:
        return 0.0
    f3 = (1.0 - nu3 / 3.0) / (1.0 - 1.0 / 3.0) ** 3
    val = 4.0 * f3 * C3
    ps = set(odd_prime_factors(d1)) | set(odd_prime_factors(d2)) | set(odd_prime_factors(abs(d1 - d2)))
    ps.discard(3)
    for p in ps:
        rr1, rr2 = d1 % p, d2 % p
        nu = len({0, rr1, rr2})
        generic = (1.0 - 3.0 / p) / (1.0 - 1.0 / p) ** 3
        actual = (1.0 - nu / p) / (1.0 - 1.0 / p) ** 3
        val *= actual / generic
    return val

print("accumulating S3 and S3pair over even pairs to HMAX ...", flush=True)
evens = np.arange(2, HMAX + 1, 2)
S2v = {int(h): S2(int(h)) for h in evens}
Hs = [250, 500, 1000, 2000, 4000]
S3cum = {H: 0.0 for H in Hs}
S3pcum = {H: 0.0 for H in Hs}
for d1 in evens:
    for d2 in evens:
        if d1 == d2:
            continue
        v = S3(int(d1), int(d2))
        vp = S2v[int(d1)] * S2v[int(d2)] * S2(abs(int(d1) - int(d2)))
        m = max(d1, d2)
        for H in Hs:
            if m <= H:
                S3cum[H] += v
                S3pcum[H] += vp
    if d1 % 400 == 0:
        print(f"  d1={d1} done", flush=True)

print("\n  H       S3(H)          S3pair(H)      S3-H^2/H     (S3-S3pair)/H")
rowsA, rowsB = [], []
for H in Hs:
    exA = (S3cum[H] - float(H) * H) / H
    exB = (S3cum[H] - S3pcum[H]) / H
    rowsA.append((H, exA)); rowsB.append((H, exB))
    print(f"  {H:5d}  {S3cum[H]:14.1f} {S3pcum[H]:14.1f} {exA:+12.4f} {exB:+12.4f}")

def fit_blog(rows):
    """fit ex(H) = b*log H + c  (the H^2 leading term already removed)."""
    xs = np.log(np.array([r[0] for r in rows], dtype=float))
    ys = np.array([r[1] for r in rows])
    return np.polyfit(xs, ys, 1)[0]

bA = fit_blog(rowsA)
bB = fit_blog(rowsB)
vA = ("PAIRWISE-DOMINATED (-1.5)" if abs(bA + 1.5) < 0.2 else
      "NULL (0)" if abs(bA) < 0.2 else f"NEITHER (b={bA:+.3f}) - flag")
vB = ("PAIRWISE-EXHAUSTED" if abs(bB) < 0.1 else f"GENUINE k=3 RIGIDITY (b={bB:+.3f})")
print(f"\n  A3: b = {bA:+.4f}  ->  {vA}")
print(f"  B3: b = {bB:+.4f}  ->  {vB}")

# ---------------- C3: psi-window skewness ----------------
X = 30_000_000
print(f"\nC3: sieving Lambda to {X} ...", flush=True)
lam = np.zeros(X + 1)
comp = np.zeros(X + 1, dtype=bool)
for p in range(2, X + 1):
    if not comp[p]:
        comp[2 * p::p] = True
        lp = np.log(p)
        pk = p
        while pk <= X:
            lam[pk] = lp
            pk *= p
cum = np.cumsum(lam)
del lam, comp
print("  H        skew      3*std_gauss   verdict")
for H in (100, 300, 1000, 3000, 10000):
    xs_grid = np.arange(X // 2, X - H, H)
    d = cum[xs_grid + H] - cum[xs_grid] - H
    d = d - np.mean(d)
    m2 = float(np.mean(d * d)); m3 = float(np.mean(d * d * d))
    skew = m3 / m2 ** 1.5
    band = 3 * np.sqrt(15.0 / len(xs_grid))
    print(f"  {H:8d}  {skew:+8.4f}   {band:8.4f}      {'EXCEEDS' if abs(skew) > band else 'gaussian'}",
          flush=True)
