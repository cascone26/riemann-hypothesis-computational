#!/usr/bin/env python3
"""§832: the capacity-matched KNIFE-EDGE ladder — does the certified below-capacity Weil margin m(4W)
stay bounded below as the function space GROWS (4W→∞ with N tracking frame capacity), or does it thin
toward 0? This is the B16/§799/§811 load-bearing question ("keep m(W) alive while dim K_W ~ e^{4W}").

Built on the §830-validated instrument: below-capacity Fejér spans certify cleanly (radii ~1e-6 ≪
margin ~0.1). Here N = round(f·N_cap), N_cap = 15·4W/π, kept safely below capacity (well-conditioned).
Certify the Weil pencil-min (Q=pole+arch−prime vs G) at each rung; report m(4W) and its trend.
Reuses the §788/§823/§830 zero-free arb entry functions directly (no module side-effect)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
ctx.prec = 250

def certmargin(fourW, CEN, rc=400):
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
    for _ in range(50):
        mid = (lo+hi)/2
        if ok(arb(mid)): lo = mid
        else: hi = mid
    return lo, maxrad

print("=" * 92)
print("§832 CAPACITY-MATCHED KNIFE-EDGE LADDER — does the certified margin m(4W) survive as the space grows?")
print("N ≈ 0.78·N_cap (N_cap=15·4W/π), window [14,44]; margin = certified Weil pencil-min (Handle if >0)")
print("=" * 92)
print(f"{'4W':>6} {'N_cap':>7} {'N':>4} {'N/Ncap':>7} {'nprime':>7} {'cert margin m(4W)':>19} {'maxrad':>9}  verdict")
rungs = [1.30, 1.60, 1.90, 2.20, 2.50]
results = []
LOGN = np.log(np.arange(2, 300))
for fourW in rungs:
    Ncap = 15 * fourW / np.pi
    N = max(4, round(0.78 * Ncap))
    CEN = list(np.linspace(14, 44, N))
    npr = int((LOGN < fourW).sum())
    m, mr = certmargin(fourW, CEN)
    if m is not None:
        results.append((fourW, m)); val = f"{m:.6f}"; verdict = "CERTIFIED (Handle)"
    else:
        val = "—"; verdict = "radii-swallowed"
    print(f"{fourW:>6.2f} {Ncap:>7.2f} {N:>4} {N/Ncap:>7.2f} {npr:>7} {val:>19} {mr:>9.1e}  {verdict}", flush=True)

print("=" * 92)
if len(results) >= 2:
    ms = [m for _, m in results]
    trend = "BOUNDED BELOW (knife-edge survives)" if min(ms) > 0.5*max(ms) else ("THINNING toward 0" if ms[-1] < ms[0] else "rising")
    print(f"certified margins m(4W): {[round(m,4) for _,m in results]}")
    print(f"min={min(ms):.4f}  max={max(ms):.4f}  first={ms[0]:.4f}  last={ms[-1]:.4f}  →  {trend}")
    print("If bounded below across the ladder: strong evidence the capacity-matched Weil margin does NOT")
    print("thin as the space grows — the knife-edge holds on the certified below-capacity envelope.")
print("=" * 92)
