#!/usr/bin/env python3
"""§874: push the HARD low-window corner. Fine arch integration (seg=0.5, prec 380) to shrink radii below the
razor-thin corner margin, certifying [14,44] past the 2.5 frontier."""
import sys,os,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from flint import arb,acb,ctx
import arb_fejer_cert2 as ac, arb_zerofree_cert3 as zc
PREC=380
def arch_fine(cj,ck,RC,seg=0.5):
    def integ(r,_): return ac.fej_acb(acb(cj),r)*ac.fej_acb(acb(ck),r)*ac.omega_acb(r)
    ctx.prec=420; I=acb(0); a=arb(-RC); RCa=arb(RC); s=arb(seg)
    while a<RCa:
        b=a+s
        if b>RCa: b=RCa
        I=I+acb.integral(integ,acb(a),acb(b)); a=b
    ctx.prec=PREC
    val=I.real/(2*arb.pi()); W=ac.W; B1=ac.BAND1
    tail=(arb(2)/arb.pi())*(arb(RC).log()+2)/(W**4*2*(arb(RC)-B1)**3)*2
    return val+arb(0).union(tail.union(-tail))
def certQ(fourW,CEN,rc):
    ctx.prec=PREC; ac.W=arb(fourW)/4; ac.TW=2*ac.W; ac.c_tri=arb.pi()/ac.W
    cs=[arb(c) for c in CEN]; n=len(cs); Q=[[None]*n for _ in range(n)]; Gm=[[None]*n for _ in range(n)]; mr=0.0
    for j in range(n):
        for k in range(j,n):
            Q[j][k]=Q[k][j]=ac.pole_entry(cs[j],cs[k])+arch_fine(cs[j],cs[k],rc)-ac.prime_entry(cs[j],cs[k])
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
    lo,hi=0.0,0.5
    for _ in range(46):
        mid=(lo+hi)/2
        if ok(arb(mid)): lo=mid
        else: hi=mid
    return lo,mr
a,b,fourW,rc=14,44,2.75,180
CEN=[float(x) for x in np.linspace(a,b,8)]
m,mr=certQ(fourW,CEN,rc)
print(f"[{a},{b}] 4W={fourW} N=8 rc={rc} (fine seg=0.5): cert Q>= {m if m is None else round(m,8)}  maxrad {mr:.2e}  {'CERTIFIED' if m else 'no (margin<radii)'}")
