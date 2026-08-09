#!/usr/bin/env python3
"""§849: extend the certified positivity set to HIGHER zero-heights. §830 certified Weil positivity on
below-capacity spans at low windows; §848 certified the arch floor GROWS with height. Combining: at higher
windows the bigger floor should make the FULL Weil margin (arch+pole−prime) certify even more robustly.
Certify Q pencil-min > 0 (full Weil form, zero-free/unconditional as in §788) at [60,90] and [120,150],
below capacity — a new certified Handle: Weil positivity for band-limited test functions at higher heights.
Reuses the §788/§823 arb entry functions directly (no module side-effect)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
ctx.prec = 220

def certQ(fourW, CEN, rc):
    ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
    cs = [arb(c) for c in CEN]; n = len(cs)
    Q = [[None]*n for _ in range(n)]; G = [[None]*n for _ in range(n)]; maxrad = 0.0
    for j in range(n):
        for k in range(j, n):
            Q[j][k] = Q[k][j] = ac.pole_entry(cs[j], cs[k]) + zc.arch_entry_v3(cs[j], cs[k], arb(rc)) - ac.prime_entry(cs[j], cs[k])
            G[j][k] = G[k][j] = 2*arb.pi()*ac.gpair_closed_arb(cs[j], cs[k], arb(0))
            maxrad = max(maxrad, float(Q[j][k].rad()))
    def ok(lam):
        d = [None]*n; L = [[arb(0)]*n for _ in range(n)]
        for i in range(n):
            sd = Q[i][i] - lam*G[i][i]
            for kk in range(i): sd -= L[i][kk]*L[i][kk]*d[kk]
            if not (sd > 0): return False
            d[i] = sd
            for j2 in range(i+1, n):
                v = Q[j2][i] - lam*G[j2][i]
                for kk in range(i): v -= L[j2][kk]*L[i][kk]*d[kk]
                L[j2][i] = v/d[i]
        return True
    if not ok(arb(0)): return None, maxrad
    lo, hi = 0.0, 1.0
    for _ in range(48):
        mid = (lo+hi)/2
        if ok(arb(mid)): lo = mid
        else: hi = mid
    return lo, maxrad

print("="*92)
print("§849 CERTIFY Weil positivity at HIGHER zero-heights (full Q=pole+arch−prime, zero-free/unconditional)")
print("="*92)
print(f"{'window':>14} {'height':>7} {'4W':>4} {'N':>3} {'cert Q pencil-min ≥':>19} {'maxrad':>9}  verdict")
for (a, b, fourW, rc) in [(60,90,2.0,300),(120,150,2.0,360),(60,90,2.5,300),(120,150,2.5,360)]:
    L = b-a; W = fourW/4; mid = (a+b)/2
    N = max(4, round(0.6*L*2*W/np.pi))
    CEN = [float(x) for x in np.linspace(a, b, N)]
    m, mr = certQ(fourW, CEN, rc)
    if m is not None:
        print(f"{f'[{a},{b}]':>14} {mid:>7.0f} {fourW:>4.1f} {N:>3} {m:>19.6f} {mr:>9.1e}  CERTIFIED (Handle)", flush=True)
    else:
        print(f"{f'[{a},{b}]':>14} {mid:>7.0f} {fourW:>4.1f} {N:>3} {'NOT PD':>19} {mr:>9.1e}  not certified", flush=True)
print("="*92)
print("Certified Q>0 = Weil positivity for band-limited test functions at these higher zero-heights,")
print("unconditional (zero-free route, no verified-zero list). Extends §808/§830 certified set to higher t.")
