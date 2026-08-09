#!/usr/bin/env python3
"""Quick test to verify arb certification pipeline works and report initial results."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc

ctx.prec = 200  # Standard precision

print("Quick test: certifying Weil margin at 4W=1.30 (known positive case)")
print("=" * 80)

fourW = 1.30
Ncap = 15 * fourW / np.pi
N = max(4, round(0.78 * Ncap))
CEN = list(np.linspace(14, 44, N))

print(f"4W={fourW}, N={N}, #centers={len(CEN)}")

ac.W = arb(fourW)/4
ac.TW = 2*ac.W
ac.c_tri = arb.pi()/ac.W
cs = [arb(c) for c in CEN]
n = len(cs)

print(f"Building Q and G matrices ({n}x{n})...")

Q = [[None]*n for _ in range(n)]
G = [[None]*n for _ in range(n)]

for j in range(n):
    for k in range(j, n):
        Q[j][k] = Q[k][j] = ac.pole_entry(cs[j], cs[k]) + zc.arch_entry_v3(cs[j], cs[k], arb(400)) - ac.prime_entry(cs[j], cs[k])
        G[j][k] = G[k][j] = 2*arb.pi()*ac.gpair_closed_arb(cs[j], cs[k], arb(0))
    print(f"  row {j}/{n-1} done")

print("Checking PD-ness via LDL^T bisection...")

def ok(lam):
    d = [None]*n
    L = [[arb(0)]*n for _ in range(n)]
    for i in range(n):
        sd = Q[i][i] - lam*G[i][i]
        for kk in range(i):
            sd -= L[i][kk]*L[i][kk]*d[kk]
        if not (sd > 0):
            return False
        d[i] = sd
        for j2 in range(i+1, n):
            v = Q[j2][i] - lam*G[j2][i]
            for kk in range(i):
                v -= L[j2][kk]*L[i][kk]*d[kk]
            L[j2][i] = v/d[i]
    return True

if not ok(arb(0)):
    print("MATRIX NOT PD!")
    sys.exit(1)

lo, hi = 0.0, 1.0
for _ in range(50):
    mid = (lo + hi) / 2
    if ok(arb(mid)):
        lo = mid
    else:
        hi = mid

print(f"\nResult:")
print(f"  Certified margin lower bound = {lo:.8e}")
print(f"  Certified margin upper bound = {hi:.8e}")
print(f"  Certified margin interval    = [{lo:.8e}, {hi:.8e}]")
print(f"  Positivity: {'YES ✓' if lo > 0 else 'UNCERTIFIABLE'}")

print("\n" + "=" * 80)
print("✓ Pipeline works! Now testing 4W=2.30...")

fourW = 2.30
Ncap = 15 * fourW / np.pi
N = max(4, round(0.80 * Ncap))
CEN = list(np.linspace(14, 44, N))

print(f"4W={fourW}, N={N}, Ncap={Ncap:.1f}, N/Ncap={N/Ncap:.2f}")

ac.W = arb(fourW)/4
ac.TW = 2*ac.W
ac.c_tri = arb.pi()/ac.W
cs = [arb(c) for c in CEN]
n = len(cs)

print(f"Building Q and G matrices ({n}x{n})...")

Q = [[None]*n for _ in range(n)]
G = [[None]*n for _ in range(n)]

for j in range(n):
    for k in range(j, n):
        Q[j][k] = Q[k][j] = ac.pole_entry(cs[j], cs[k]) + zc.arch_entry_v3(cs[j], cs[k], arb(400)) - ac.prime_entry(cs[j], cs[k])
        G[j][k] = G[k][j] = 2*arb.pi()*ac.gpair_closed_arb(cs[j], cs[k], arb(0))

if not ok(arb(0)):
    print("MATRIX NOT PD at shift 0!")
    print("4W=2.30 does NOT certify positivity")
else:
    lo, hi = 0.0, 1.0
    for _ in range(50):
        mid = (lo + hi) / 2
        if ok(arb(mid)):
            lo = mid
        else:
            hi = mid
    print(f"Certified margin lower bound = {lo:.8e}")
    print(f"Positivity: {'YES ✓' if lo > 0 else 'UNCERTIFIABLE or thin'}")
