#!/usr/bin/env python3
"""CERT CAMPAIGN (§824): sustained certification of the §814-821 structural picture. Three certified
curves, logged to cert_campaign.txt as they land:
 (A) arch_floor(W)  — certified min-eig of pencil(arch,G), the fixed archimedean budget.
 (B) frame_c(W)     — certified frame constant Σ_{γ≤800}ψ̂(γ)²≥c‖ψ‖² ABOVE t* (the provable half).
 (C) Qpos(W)        — certified Weil positivity Q_pencil_min, with instrument-horizon flag.
All arb, prec 250, fixed-segment (no hang), module state set. Runs long; append+flush per point."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
import math
OUT="cert_campaign.txt"
def setW(fourW): ac.W=arb(fourW)/4; ac.TW=2*ac.W; ac.c_tri=arb.pi()/ac.W
def pencil_min(M,G,n):
    def ok(lam):
        d=[None]*n; L=[[arb(0)]*n for _ in range(n)]
        for i in range(n):
            sd=M[i][i]-lam*G[i][i]
            for kk in range(i): sd-=L[i][kk]*L[i][kk]*d[kk]
            if not (sd>0): return False
            d[i]=sd
            for j2 in range(i+1,n):
                v=M[j2][i]-lam*G[j2][i]
                for kk in range(i): v-=L[j2][kk]*L[i][kk]*d[kk]
                L[j2][i]=v/d[i]
        return True
    if not ok(arb(0)): return None
    lo,hi=0.0,1.5
    for _ in range(44):
        mid=(lo+hi)/2
        lo,hi=(mid,hi) if ok(arb(mid)) else (lo,mid)
    return lo
def Gmat(cs,n): return [[2*arb.pi()*ac.gpair_closed_arb(cs[j],cs[k],arb(0)) for k in range(n)] for j in range(n)]
def span(lo,hi,npr): return [arb(round(lo+(hi-lo)*i/(npr-1),3)) for i in range(npr)]
def log(s): open(OUT,"a").write(s+"\n"); print(s,flush=True)

ctx.prec=250
log("# CERT CAMPAIGN §824 — arch_floor / frame_c(above t*) / Qpos, all arb-certified")
# (A) arch floor curve
log("## (A) arch_floor(W) certified:")
for fourW in [1.20,1.35,1.50,1.65,1.80,1.95,2.10,2.30]:
    setW(fourW); npr=max(4,min(6,round(30*(fourW/2)/math.pi/1.3)))
    cs=span(14,44,npr); n=len(cs)
    ARCH=[[zc.arch_entry_v3(cs[j],cs[k],arb(400)) for k in range(n)] for j in range(n)]
    m=pencil_min(ARCH,Gmat(cs,n),n)
    log(f"A {fourW:.2f} npr={npr} arch_floor {'>= %.5f'%m if m else 'NOT PD'}")
# (B) frame constant above t*
log("## (B) frame_c(W) above t* certified (zeros<=800):")
ZC=np.load('zeros500.npy'); zeros=ZC[ZC<=800.0]
def fej_arb(c,r,W):
    o=arb(0)
    for cc in (c,-c):
        y=W*(r-cc); o=o+(arb(1).union((y.sin()/y)**2) if y.contains(arb(0)) else (y.sin()/y)**2)
    return o
for fourW in [1.00,1.32,1.61,2.00]:
    setW(fourW); ts=2*math.pi*math.exp(fourW); lo=math.ceil(ts+2)
    npr=4; cs=span(lo,lo+6*(npr-1),npr); n=len(cs)  # capacity span above t*
    if float(cs[-1])>44+30: cs=span(lo,44+30,npr)
    S=[[arb(0)]*n for _ in range(n)]
    for g in zeros:
        if g>800: break
        gb=arb(float(g),1e-9); F=[fej_arb(cs[j],gb,ac.W) for j in range(n)]
        for j in range(n):
            for k in range(j,n): S[j][k]=S[j][k]+F[j]*F[k]
    for j in range(n):
        for k in range(j,n): S[k][j]=S[j][k]
    m=pencil_min([[2*S[j][k] for k in range(n)] for j in range(n)],Gmat(cs,n),n)
    log(f"B {fourW:.2f} t*={ts:.1f} span[{float(cs[0]):.0f},{float(cs[-1]):.0f}] frame_c {'>= %.5f'%m if m else 'NOT PD'}")
# (C) Q positivity across crossover, instrument-horizon flagged
log("## (C) Qpos(W) certified (Weil positivity, horizon-flagged):")
for fourW in [1.32,1.50,1.70,1.90,2.10]:
    setW(fourW); npr=max(4,min(6,round(30*(fourW/2)/math.pi/1.3)))
    cs=span(14,44,npr); n=len(cs)
    Q=[[ac.pole_entry(cs[j],cs[k])+zc.arch_entry_v3(cs[j],cs[k],arb(400))-ac.prime_entry(cs[j],cs[k]) for k in range(n)] for j in range(n)]
    mr=max(float(Q[j][j].rad()) for j in range(n))
    m=pencil_min(Q,Gmat(cs,n),n)
    log(f"C {fourW:.2f} npr={npr} maxrad={mr:.1e} Qpos {'>= %.6f CERTIFIED'%m if m else 'INSTRUMENT-UNDECIDABLE'}")
log("# CAMPAIGN COMPLETE")
