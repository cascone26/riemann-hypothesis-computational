#!/usr/bin/env python3
"""§831: the precision-frontier probe — is the critical-sampling margin collapse a NUMERICAL artifact
(higher arb precision rescues it) or a GENUINE degeneracy (margin truly →0, no precision helps)?

This distinguishes the two readings §829 left open. Fix 4W=1.60 (N_cap=30·2W/π≈7.6). For N below,
at, and past capacity, crank ctx.prec (100→400 bits) and at each precision certify the Weil pencil-min
(Q=pole+arch−prime vs G) via interval LDLᵀ bisection. Record: does the certified margin STABILIZE at a
positive value as precision rises (→ real positive margin, just needed resolution), or does it stay
radii-swallowed / →0 no matter the precision (→ genuine critical-sampling degeneracy = the §823 razor
is real, not instrumental)?

Reuses the §788/§823 zero-free arb machinery via cert_subtstar's entry functions (module state set)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc

def certQ_prec(fourW, CEN, prec, rc=400):
    """Certified pencil-min of (Q,G) at a given arb precision. Returns (margin_or_None, maxrad, minGrad)."""
    ctx.prec = prec
    ac.W = arb(fourW) / 4; ac.TW = 2 * ac.W; ac.c_tri = arb.pi() / ac.W
    cs = [arb(c) for c in CEN]; n = len(cs)
    Q = [[None]*n for _ in range(n)]; G = [[None]*n for _ in range(n)]; maxrad = 0.0
    for j in range(n):
        for k in range(j, n):
            Q[j][k] = Q[k][j] = ac.pole_entry(cs[j], cs[k]) + zc.arch_entry_v3(cs[j], cs[k], arb(rc)) - ac.prime_entry(cs[j], cs[k])
            G[j][k] = G[k][j] = 2*arb.pi()*ac.gpair_closed_arb(cs[j], cs[k], arb(0))
            maxrad = max(maxrad, float(Q[j][k].rad()))
    # min diagonal of G as a crude conditioning proxy (radius-aware)
    minGrad = min(float(G[i][i].rad()) for i in range(n))
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
    if not ok(arb(0)): return None, maxrad, minGrad
    lo, hi = 0.0, 1.0
    for _ in range(50):
        mid = (lo+hi)/2
        if ok(arb(mid)): lo = mid
        else: hi = mid
    return lo, maxrad, minGrad

fourW = 1.60
Ncap = 30 * (2*(fourW/4)) / np.pi
print("=" * 96)
print(f"§831 PRECISION-FRONTIER PROBE at 4W={fourW} (N_cap=30·2W/π≈{Ncap:.2f}); t*={2*np.pi*np.exp(fourW):.1f}")
print("Q: does higher precision RESCUE the margin (numerical) or NOT (genuine critical-sampling)?")
print("=" * 96)
for N in [6, 7, 8, 9]:
    CEN = list(np.linspace(14, 44, N))
    reg = "below cap" if N < Ncap-0.5 else ("~at cap" if N < Ncap+1.5 else "OVER cap")
    print(f"\nN={N} ({reg}), centers {[round(c,1) for c in CEN]}:")
    print(f"  {'prec':>6} {'cert margin':>16} {'maxrad':>10}")
    prev = None
    for prec in [100, 200, 300, 400]:
        m, mr, mg = certQ_prec(fourW, CEN, prec)
        val = f"{m:.3e}" if m is not None else "radii-swallowed"
        print(f"  {prec:>6} {val:>16} {mr:>10.1e}", flush=True)
print("\n" + "=" * 96)
print("READ: margin STABILIZES to a positive plateau as prec rises → real positive margin (resolution-limited).")
print("      margin stays radii-swallowed / keeps shrinking with prec → genuine degeneracy (razor is real).")
print("=" * 96)
