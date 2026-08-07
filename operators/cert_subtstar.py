#!/usr/bin/env python3
"""§823: does certification REACH the sub-t* corner? Certify Q=arch+pole−prime PD (Weil positivity)
on a DENSE capacity span at 4W=1.90 (deep past crossover t*≈42; float margin ~0.0008 razor-thin).
Certifies → certified sub-t* positivity. Radii swallow it → instrument-undecidable (Gödel §806).
Uses the §788 zero-free machinery (module state set)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
ctx.prec=250
def certQ(fourW, CEN, rc=400):
    ac.W=arb(fourW)/4; ac.TW=2*ac.W; ac.c_tri=arb.pi()/ac.W
    cs=[arb(c) for c in CEN]; n=len(cs)
    Q=[[None]*n for _ in range(n)]; G=[[None]*n for _ in range(n)]; maxrad=0.0
    for j in range(n):
        for k in range(j,n):
            Q[j][k]=Q[k][j]=ac.pole_entry(cs[j],cs[k])+zc.arch_entry_v3(cs[j],cs[k],arb(rc))-ac.prime_entry(cs[j],cs[k])
            G[j][k]=G[k][j]=2*arb.pi()*ac.gpair_closed_arb(cs[j],cs[k],arb(0))
            maxrad=max(maxrad,float(Q[j][k].rad()))
    def ok(lam):
        d=[None]*n; L=[[arb(0)]*n for _ in range(n)]
        for i in range(n):
            sd=Q[i][i]-lam*G[i][i]
            for kk in range(i): sd-=L[i][kk]*L[i][kk]*d[kk]
            if not (sd>0): return False
            d[i]=sd
            for j2 in range(i+1,n):
                v=Q[j2][i]-lam*G[j2][i]
                for kk in range(i): v-=L[j2][kk]*L[i][kk]*d[kk]
                L[j2][i]=v/d[i]
        return True
    pd=ok(arb(0))
    if not pd: return None, maxrad
    lo,hi=0.0,1.0
    for _ in range(45):
        mid=(lo+hi)/2
        if ok(arb(mid)): lo=mid
        else: hi=mid
    return lo, maxrad
# sub-t* dense span at 4W=1.90 (t*~42, whole window [14,44] near/below t*)
for fourW,CEN in [(1.90,[16,22,28,34,40]),(1.90,[15,20,25,30,35,40]),(2.00,[16,22,28,34,40])]:
    m,mr=certQ(fourW,CEN)
    tag=f"CERTIFIED Q_pencil_min ≥ {m:.6f} (sub-t* Weil positivity CERTIFIED)" if m else f"NOT certified (radii {mr:.1e} vs razor margin — INSTRUMENT-UNDECIDABLE)"
    print(f"4W={fourW} span {CEN} (maxrad {mr:.1e}): {tag}", flush=True)
