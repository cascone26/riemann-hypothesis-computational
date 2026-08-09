#!/usr/bin/env python3
import sys,os,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from flint import arb,ctx
import arb_fejer_cert2 as ac, arb_zerofree_cert3 as zc
ctx.prec=170
def certQ(fourW,CEN,rc):
    ac.W=arb(fourW)/4; ac.TW=2*ac.W; ac.c_tri=arb.pi()/ac.W
    cs=[arb(c) for c in CEN]; n=len(cs); Q=[[None]*n for _ in range(n)]; Gm=[[None]*n for _ in range(n)]; mr=0.0
    for j in range(n):
        for k in range(j,n):
            Q[j][k]=Q[k][j]=ac.pole_entry(cs[j],cs[k])+zc.arch_entry_v3(cs[j],cs[k],arb(rc))-ac.prime_entry(cs[j],cs[k])
            Gm[j][k]=Gm[k][j]=2*arb.pi()*ac.gpair_closed_arb(cs[j],cs[k],arb(0))
            mr=max(mr,float(Q[j][k].rad()))
    def ok(lam):
        d=[None]*n; L=[[arb(0)]*n for _ in range(n)]
        for i in range(n):
            sd=Q[i][i]-lam*Gm[i][i]
            for kk in range(i): sd-=L[i][kk]*L[i][kk]*d[kk]
            if not (sd>0): return False
            d[i]=sd
            for j2 in range(i+1,n):
                v=Q[j2][i]-lam*Gm[j2][i]
                for kk in range(i): v-=L[j2][kk]*L[i][kk]*d[kk]
                L[j2][i]=v/d[i]
        return True
    if not ok(arb(0)): return None,mr
    lo,hi=0.0,1.0
    for _ in range(44):
        mid=(lo+hi)/2
        if ok(arb(mid)): lo=mid
        else: hi=mid
    return lo,mr
a,b,fourW,rc=120,150,4.0,300
CEN=[float(x) for x in np.linspace(a,b,4)]
m,mr=certQ(fourW,CEN,rc)
print(f"[{a},{b}] 4W={fourW} N=4 rc={rc}: cert Q>= {m if m is None else round(m,6)}  maxrad {mr:.1e}  {'CERTIFIED' if m else 'no'}")
