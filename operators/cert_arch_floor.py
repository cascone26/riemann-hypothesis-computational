#!/usr/bin/env python3
"""§816-CERT: arb-certify the archimedean-floor W-invariance. Certified min-eig of the arch pencil
(arch, G) on a capacity span, at several W — turns the §816 float observation (arch floor ~0.153
W-invariant) into a certified structural HANDLE. Uses arch_entry_v3 (arb, fixed-segment, no hang)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
ctx.prec=200
def arch_floor_cert(fourW, CEN):
    ac.W=arb(fourW)/4; ac.TW=2*ac.W; ac.c_tri=arb.pi()/ac.W
    cs=[arb(c) for c in CEN]; n=len(cs)
    ARCH=[[None]*n for _ in range(n)]; G=[[None]*n for _ in range(n)]
    for j in range(n):
        for k in range(j,n):
            ARCH[j][k]=ARCH[k][j]=zc.arch_entry_v3(cs[j],cs[k],arb(400))
            G[j][k]=G[k][j]=2*arb.pi()*ac.gpair_closed_arb(cs[j],cs[k],arb(0))
    # certified generalized min-eig of (ARCH, G): bisect λ, interval-LDLᵀ of (ARCH-λG)
    def ok(lam):
        d=[None]*n; L=[[arb(0)]*n for _ in range(n)]
        for i in range(n):
            sd=ARCH[i][i]-lam*G[i][i]
            for kk in range(i): sd-=L[i][kk]*L[i][kk]*d[kk]
            if not (sd>0): return False
            d[i]=sd
            for j2 in range(i+1,n):
                v=ARCH[j2][i]-lam*G[j2][i]
                for kk in range(i): v-=L[j2][kk]*L[i][kk]*d[kk]
                L[j2][i]=v/d[i]
        return True
    if not ok(arb(0)): return None
    lo,hi=0.0,1.0
    for _ in range(45):
        mid=(lo+hi)/2
        if ok(arb(mid)): lo=mid
        else: hi=mid
    return lo
print("4W   span            certified arch_floor (min-eig of pencil (arch,G))")
for fourW,CEN in [(1.32,[20,26,32,38]),(1.60,[18,24,30,36,42]),(1.90,[16,22,28,34,40])]:
    af=arch_floor_cert(fourW,CEN)
    print(f"{fourW}  {CEN}  {'CERTIFIED arch_floor ≥ %.5f'%af if af else 'not PD'}", flush=True)
