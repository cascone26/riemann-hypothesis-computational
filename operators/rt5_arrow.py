#!/usr/bin/env python3
"""§352: THE k=2 ARROW EXPERIMENTS — pre-registered in LEDGER §351 BEFORE this run.

A: slope of [sum_{h<=H} S(h) - H] vs log H.       NULL |slope|<0.05 | ARROW -0.50+-0.05
B: slope of [Var(psi(x+H)-psi(x)-H)/H - log xbar] vs log H.  NULL 0+-0.1 | ARROW -1.0+-0.15
C: our DT entry vs Poisson / Hardy-Littlewood covariance predictions.
   PASS: actual^2 <= 3x HL and HL <= Poisson/3.
S(h) = 2 C2 prod_{p|h, p>2} (p-1)/(p-2) for even h, 0 for odd. All zero-free.
"""
import numpy as np

# ---------------- Part A: singular series sum ----------------
HMAX = 2_000_000
print("Part A: sieving smallest-prime-factor to", HMAX, flush=True)
spf = np.zeros(HMAX + 1, dtype=np.int32)
for p in range(2, int(HMAX ** 0.5) + 1):
    if spf[p] == 0:
        sl = spf[p * p::p]
        sl[sl == 0] = p
        spf[p * p::p] = sl
# C2 from primes up to 1e6
PR = np.array([p for p in range(3, 1_000_000) if (spf[p] == 0 if p <= HMAX else True) and p % 2 == 1
               and (spf[p] == 0)][:0])  # placeholder avoided; compute below properly
primes_mask = spf[: 1_000_000 + 1] == 0
primes_mask[:2] = False
primes = np.nonzero(primes_mask)[0]
odd_primes = primes[primes > 2].astype(np.float64)
C2 = float(np.prod(1.0 - 1.0 / (odd_primes[:200000] - 1.0) ** 2))
print(f"  C2 (twin-prime const, computed) = {C2:.10f}")

S = np.zeros(HMAX + 1)
for h in range(2, HMAX + 1, 2):
    m = h
    prod = 1.0
    # strip 2s
    while m % 2 == 0:
        m //= 2
    while m > 1:
        p = spf[m] if spf[m] != 0 else m
        prod *= (p - 1.0) / (p - 2.0)
        while m % p == 0:
            m //= p
    S[h] = 2.0 * C2 * prod
CS = np.cumsum(S)
print("  H        sum S(h)-H      [per-decade slope vs log H]")
Hs = [1000, 3000, 10000, 30000, 100000, 300000, 1000000, 2000000]
prev = None
slopes = []
for H in Hs:
    d = CS[H] - H
    line = f"  {H:8d}  {d:+12.4f}"
    if prev is not None:
        sl = (d - prev[1]) / (np.log(H) - np.log(prev[0]))
        slopes.append(sl)
        line += f"   {sl:+.4f}"
    prev = (H, d)
    print(line, flush=True)
# least-squares over the full range
xs = np.log(np.array(Hs, dtype=float))
ys = np.array([CS[H] - H for H in Hs])
A_ls = np.polyfit(xs, ys, 1)
slopeA = A_ls[0]
verdictA = ("ARROW (-0.5)" if abs(slopeA + 0.5) < 0.05 else
            "NULL (0)" if abs(slopeA) < 0.05 else "NEITHER - flag")
print(f"  Part A least-squares slope = {slopeA:+.4f}  ->  {verdictA}")

# ---------------- Part B: psi short-interval variance ----------------
X = 30_000_000
print(f"\nPart B: sieving Lambda to {X} ...", flush=True)
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
print("  H        Var/H - log(xbar)   [slope vs log H]")
res = []
for H in (100, 300, 1000, 3000, 10000, 30000, 100000):
    xs_grid = np.arange(X // 2, X - H, H)          # non-overlapping
    d = cum[xs_grid + H] - cum[xs_grid] - H
    var = float(np.mean(d * d))
    xbar = float(np.mean(xs_grid)) + H / 2
    val = var / H - np.log(xbar)
    res.append((H, val))
    line = f"  {H:8d}  {val:+12.4f}"
    if len(res) > 1:
        sl = (res[-1][1] - res[-2][1]) / (np.log(res[-1][0]) - np.log(res[-2][0]))
        line += f"   {sl:+.4f}"
    print(line, flush=True)
xsb = np.log(np.array([r[0] for r in res], dtype=float))
ysb = np.array([r[1] for r in res])
slopeB = np.polyfit(xsb, ysb, 1)[0]
verdictB = ("ARROW (-1)" if abs(slopeB + 1.0) < 0.15 else
            "NULL (0)" if abs(slopeB) < 0.1 else "NEITHER - flag")
print(f"  Part B least-squares slope = {slopeB:+.4f}  ->  {verdictB}")

# ---------------- Part C: the bridge to our DT ----------------
print("\nPart C: model-predicted E[DT_entry^2] for the [14,44] w=1.5 edge window, n0=100")
import l2_edge_sweep as L
w = 1.5
NC = 5000
n0 = 100
cs = L.probes(L.BAND, 4.4)      # a certified-region spacing (§346 reach at n0=100, w=1.5)
cj = ck = cs[0]                  # edge probe diagonal entry
ns = np.arange(n0, NC)
f_n = 2.0 * L.g_pair(cj, ck, w, np.log(ns)) / np.sqrt(ns)
# actual DT entry over [n0, NC): sum Lambda(n) f - int f dx (dx via fine grid)
lamC = np.zeros(NC)
compC = np.zeros(NC, dtype=bool)
for p in range(2, NC):
    if not compC[p]:
        compC[2 * p::p] = True
        lp = np.log(p)
        pk = p
        while pk < NC:
            lamC[pk] = lp
            pk *= p
actual = float(np.sum(lamC[n0:NC] * f_n) )
xg = np.arange(n0, NC, 0.01)
smooth = float(np.sum(2.0 * L.g_pair(cj, ck, w, np.log(xg)) / np.sqrt(xg)) * 0.01)
DTa = actual - smooth
# Poisson: sum f^2 * log n   (Var[Lambda(n)] ~ log n under Cramer)
poiss = float(np.sum(f_n * f_n * np.log(ns)))
# HL: add offdiag (S(n-m)-1); S from Part A array (covers |n-m| < NC)
Sd = S[:NC].copy()               # S(d) for d>=2 even; S(odd)=0; define S(0) unused, S(1)=0
hl = poiss
# correlation sum: sum_{d>=1} (S(d)-1) * 2 * sum_n f(n) f(n+d)
for dd in range(1, NC - n0):
    c = Sd[dd] - 1.0
    if abs(c) < 1e-12 and dd % 2 == 0:
        continue
    ov = float(np.sum(f_n[:-dd] * f_n[dd:])) if dd < len(f_n) else 0.0
    hl += 2.0 * c * ov
print(f"  actual DT entry        = {DTa:+.4e}   actual^2 = {DTa*DTa:.4e}")
print(f"  Poisson prediction     = {poiss:.4e}")
print(f"  Hardy-Littlewood pred  = {hl:.4e}   (ratio HL/Poisson = {hl/poiss:.3f})")
passC = (DTa * DTa <= 3 * hl) and (hl <= poiss / 3)
print(f"  pre-registered PASS requires: actual^2 <= 3xHL AND HL <= Poisson/3  ->  "
      f"{'PASS' if passC else 'FAIL'}")
print(f"    (actual^2 <= 3HL: {DTa*DTa <= 3*hl};  HL <= Poisson/3: {hl <= poiss/3})")
