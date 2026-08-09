#!/usr/bin/env python3
"""§848: arb-CERTIFY the floor law (the unconditional BUDGET half of the arch-vs-prime race, §843).
floor(window) = λmin(arch+pole, G) is a POSITIVE-DEFINITE quantity (no prime term) ~0.13–0.66, far above
the arb radii — so it certifies cleanly. Certify floor ≥ (certified lower bound) at sample heights and
compare to (1/2π)ln(t/2π). This is unconditional (no RH, no zeros) — a certified Handle extending §808/§830.
Reuses the §788/§823 arb arch+pole entry functions (module state set)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
ctx.prec = 220

def cert_floor(fourW, CEN, rc):
    ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
    cs = [arb(c) for c in CEN]; n = len(cs)
    AP = [[None]*n for _ in range(n)]; G = [[None]*n for _ in range(n)]; maxrad = 0.0
    for j in range(n):
        for k in range(j, n):
            AP[j][k] = AP[k][j] = ac.pole_entry(cs[j], cs[k]) + zc.arch_entry_v3(cs[j], cs[k], arb(rc))
            G[j][k] = G[k][j] = 2*arb.pi()*ac.gpair_closed_arb(cs[j], cs[k], arb(0))
            maxrad = max(maxrad, float(AP[j][k].rad()))
    def ok(lam):
        d = [None]*n; L = [[arb(0)]*n for _ in range(n)]
        for i in range(n):
            sd = AP[i][i] - lam*G[i][i]
            for kk in range(i): sd -= L[i][kk]*L[i][kk]*d[kk]
            if not (sd > 0): return False
            d[i] = sd
            for j2 in range(i+1, n):
                v = AP[j2][i] - lam*G[j2][i]
                for kk in range(i): v -= L[j2][kk]*L[i][kk]*d[kk]
                L[j2][i] = v/d[i]
        return True
    if not ok(arb(0)): return None, maxrad
    lo, hi = 0.0, 1.5
    for _ in range(48):
        mid = (lo+hi)/2
        if ok(arb(mid)): lo = mid
        else: hi = mid
    return lo, maxrad

print("="*96)
print("§848 arb-CERTIFY the floor law (unconditional budget half): floor=λmin(arch+pole,G) ≥ ?  vs (1/2π)ln(t/2π)")
print("="*96)
inv2pi = 1/(2*np.pi)
print(f"{'window':>14} {'height':>7} {'4W':>4} {'N':>3} {'cert floor ≥':>13} {'maxrad':>9} {'(1/2π)ln(t/2π)':>14} {'ratio':>7}")
for (a, b, fourW, rc) in [(14,44,2.0,300),(60,90,2.0,300),(120,150,2.0,360)]:
    L = b-a; W = fourW/4; mid = (a+b)/2
    N = max(4, round(0.6*L*2*W/np.pi))
    CEN = [float(x) for x in np.linspace(a, b, N)]
    m, mr = cert_floor(fourW, CEN, rc)
    pred = inv2pi*np.log(mid/(2*np.pi))
    if m is not None:
        print(f"{f'[{a},{b}]':>14} {mid:>7.0f} {fourW:>4.1f} {N:>3} {m:>13.6f} {mr:>9.1e} {pred:>14.6f} {m/pred:>7.3f}", flush=True)
    else:
        print(f"{f'[{a},{b}]':>14} {mid:>7.0f} {fourW:>4.1f} {N:>3} {'NOT PD':>13} {mr:>9.1e} {pred:>14.6f}", flush=True)
print("="*96)
print("Certified floor ≥ bound, UNCONDITIONAL (no RH, no zeros). If ratio→~1 the (1/2π)ln(t/2π) law is confirmed")
print("in arb; this is the budget-side Handle extending §808/§830.")
