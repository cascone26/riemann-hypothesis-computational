#!/usr/bin/env python3
"""THE m(W) CURVE (§807) — certified pencil λ_min(Q,G) at capacity-matched spans across support
widths 4W, on the fixed window [14,44]. Answers the B16 knife-edge / Grothendieck crux: does the
p-block margin m(W) stay bounded below as W grows, or collapse? Uses the §802A machinery (zero-free
Q, L² pencil G, bisected interval-LDLᵀ). Per 4W: pick probe count ~capacity, back off until G is PD
(conditioned), certify m(W). ALL module state set."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
ctx.prec = 220

def mW(fourW, nprobes):
    ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
    CEN = [round(14+30*i/(nprobes-1),3) for i in range(nprobes)]
    cs = [arb(c) for c in CEN]; n = len(cs)
    Q = [[None]*n for _ in range(n)]
    for j in range(n):
        for k in range(j,n):
            Q[j][k]=Q[k][j]=ac.pole_entry(cs[j],cs[k])+zc.arch_entry_v3(cs[j],cs[k],arb(400))-ac.prime_entry(cs[j],cs[k])
    G = [[2*arb.pi()*ac.gpair_closed_arb(cs[j],cs[k],arb(0)) for k in range(n)] for j in range(n)]
    # check G PD (conditioning) via its LDLᵀ pivots
    dG=[None]*n; LG=[[arb(0)]*n for _ in range(n)]
    for i in range(n):
        sd=G[i][i]
        for k in range(i): sd-=LG[i][k]*LG[i][k]*dG[k]
        if not (sd>0): return None, None, "G_notPD"
        dG[i]=sd
        for j2 in range(i+1,n):
            v=G[j2][i]
            for k in range(i): v-=LG[j2][k]*LG[i][k]*dG[k]
            LG[j2][i]=v/dG[i]
    gpiv=min(float(dG[i].mid()) for i in range(n))
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
    if not ok(arb(0)): return None, gpiv, "Q_notPD"
    lo,hi=0.0,2.0
    for _ in range(40):
        mid=(lo+hi)/2
        if ok(arb(mid)): lo=mid
        else: hi=mid
    return lo, gpiv, "ok"

import math
print("4W  | nprobes | Gpiv | certified m(W)=λ_min(Q,G)", flush=True)
for fourW in [1.32, 1.60, 2.00, 2.50, 3.00]:
    nstart = max(6, round(6*fourW/1.32))
    for nprobes in range(nstart, 4, -1):   # back off until G PD
        m, gpiv, status = mW(fourW, nprobes)
        if status == "ok":
            print(f"{fourW:.2f} | {nprobes:2d} | {gpiv:.2e} | m(W) = {m:.6f}", flush=True)
            break
        elif status == "Q_notPD":
            print(f"{fourW:.2f} | {nprobes:2d} | {gpiv:.2e} | Q NOT PD on capacity span — margin gone here", flush=True)
            break
        # G_notPD: too many probes, back off
    else:
        print(f"{fourW:.2f} | no conditioned span found", flush=True)
