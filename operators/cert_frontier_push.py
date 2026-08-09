#!/usr/bin/env python3
"""§870: PUSH the rigorous arb-certified positivity frontier. Workflow got arb to 4W~2.5; float margins are
COMFORTABLE (>0.08) at higher windows even at 4W=4. arb-CERTIFY Weil positivity (Q=pole+arch-prime PSD, zero-free
/unconditional) at 4W=3,4 on windows with comfortable margin -> extend the rigorous frontier far past 2.5."""
import sys,os,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from flint import arb,ctx
import arb_fejer_cert2 as ac, arb_zerofree_cert3 as zc
ctx.prec=220
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
    for _ in range(46):
        mid=(lo+hi)/2
        if ok(arb(mid)): lo=mid
        else: hi=mid
    return lo,mr
print("="*84); print("§870 PUSH rigorous frontier: arb-certify Weil positivity at 4W=3,4 (comfortable-margin windows)"); print("="*84)
print(f"{'window':>12} {'4W':>4} {'N':>3} {'cert Q>=':>12} {'maxrad':>9}  verdict")
for (a,b,fourW,rc) in [(50,80,3.0,340),(120,150,4.0,420),(300,330,4.0,620)]:
    L=b-a; W=fourW/4; N=max(4,round(0.6*L*2*W/np.pi)); CEN=[float(x) for x in np.linspace(a,b,N)]
    m,mr=certQ(fourW,CEN,rc)
    v="CERTIFIED (Handle)" if m is not None else "not certified"
    print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {N:>3} {(f'{m:.6f}' if m is not None else '--'):>12} {mr:>9.1e}  {v}",flush=True)
print("="*84)
