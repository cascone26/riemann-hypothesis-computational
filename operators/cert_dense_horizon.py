#!/usr/bin/env python3
"""§826: THE INSTRUMENT-HORIZON TEST. Certify Q on the DENSE capacity span at 4W where the margin is
razor-thin (~1e-3), at escalating precision. Either certifies (razor-thin sub-t* positivity at the very
difficulty) or goes INSTRUMENT-UNDECIDABLE (Gödel §806 horizon confirmed at the open core). Dense span
= capacity-matched (K~capacity modes, the boundary shell §812)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
import math
def certQ_dense(fourW, prec):
    ctx.prec=prec
    ac.W=arb(fourW)/4; ac.TW=2*ac.W; ac.c_tri=arb.pi()/ac.W
    cap=30*(fourW/2)/math.pi
    npr=max(8,round(cap)+2)                    # DENSE, capacity-matched (over-resolves → razor margin)
    CEN=[round(14+30*i/(npr-1),3) for i in range(npr)]; cs=[arb(c) for c in CEN]; n=len(cs)
    Q=[[None]*n for _ in range(n)]; G=[[None]*n for _ in range(n)]
    for j in range(n):
        for k in range(j,n):
            Q[j][k]=Q[k][j]=ac.pole_entry(cs[j],cs[k])+zc.arch_entry_v3(cs[j],cs[k],arb(400))-ac.prime_entry(cs[j],cs[k])
            G[j][k]=G[k][j]=2*arb.pi()*ac.gpair_closed_arb(cs[j],cs[k],arb(0))
    mr=max(float(Q[j][j].rad()) for j in range(n))
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
    if not pd: return npr,mr,None
    lo,hi=0.0,0.5
    for _ in range(50):
        mid=(lo+hi)/2
        lo,hi=(mid,hi) if ok(arb(mid)) else (lo,mid)
    return npr,mr,lo
for fourW in [1.90]:
    for prec in [300,500,700]:
        npr,mr,m=certQ_dense(fourW,prec)
        tag=f"CERTIFIED Q_min ≥ {m:.2e} (razor-thin dense sub-t* positivity CERTIFIED!)" if m else f"INSTRUMENT-UNDECIDABLE (radii {mr:.1e} swallow the margin)"
        print(f"4W={fourW} DENSE npr={npr} prec={prec} maxrad={mr:.1e}: {tag}", flush=True)
        if m: break   # certified at this precision, done
