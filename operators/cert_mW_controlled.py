#!/usr/bin/env python3
"""m(W) CONTROLLED (§807b): FIXED well-separated span {20,26,32,38} (4 probes, spacing 6, well
inside [14,44], off the Nyquist edge), certified zero-free pencil λ_min(Q,G) across 4W, prec 300
+ finer arch subdivision to shrink radii (push the §807 instrument horizon back). Maximally
controlled: only 4W varies, span/count/window all fixed. Distinguishes a true m(W) trend from the
density/edge confounds of v1. (Q = arch+pole−prime = true Weil form = PSD always; λ_min>0 = margin.)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
ctx.prec = 300
CEN = [20,26,32,38]; n = len(CEN)
def mW(fourW):
    ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
    cs = [arb(c) for c in CEN]
    Q=[[None]*n for _ in range(n)]; G=[[None]*n for _ in range(n)]
    for j in range(n):
        for k in range(j,n):
            Q[j][k]=Q[k][j]=ac.pole_entry(cs[j],cs[k])+zc.arch_entry_v3(cs[j],cs[k],arb(600))-ac.prime_entry(cs[j],cs[k])
            G[j][k]=G[k][j]=2*arb.pi()*ac.gpair_closed_arb(cs[j],cs[k],arb(0))
    def ok(lam):
        d=[None]*n; L=[[arb(0)]*n for _ in range(n)]
        for i in range(n):
            sd=Q[i][i]-lam*G[i][i]
            for k in range(i): sd-=L[i][k]*L[i][k]*d[k]
            if not (sd>0): return False
            d[i]=sd
            for j2 in range(i+1,n):
                v=Q[j2][i]-lam*G[j2][i]
                for k in range(i): v-=L[j2][k]*L[i][k]*d[k]
                L[j2][i]=v/d[i]
        return True
    maxrad=max(float(Q[i][i].rad()) for i in range(n))
    if not ok(arb(0)): return None, maxrad
    lo,hi=0.0,3.0
    for _ in range(45):
        mid=(lo+hi)/2
        if ok(arb(mid)): lo=mid
        else: hi=mid
    return lo, maxrad
print("4W  | Q maxrad | certified m(W)=λ_min(Q,G) on FIXED span {20,26,32,38}", flush=True)
for fourW in [1.32,1.60,2.00,2.50,3.00]:
    m,mr=mW(fourW)
    if m is None: print(f"{fourW:.2f} | {mr:.1e} | UNDECIDABLE (radii swallow margin at prec 300)", flush=True)
    else: print(f"{fourW:.2f} | {mr:.1e} | m(W) = {m:.6f}", flush=True)
