#!/usr/bin/env python3
"""Certified λ₇: scan probe count k on [14,44]; compression = tr(G⁻¹A) via DIRECT interval LDLᵀ
solve (no orthonormalization, no eigenvector projection → no blowup). Ky-Fan: tr(G⁻¹A) over the
k-probe span ≤ Σ_top-k(true) when... NO — tr(G⁻¹A) over k-dim span = Σ of all k pencil eigenvalues.
For a k-dim span this = the k concentration values of the span; as k grows toward/past capacity it
approaches Σ_top-k(true). Need compression > tr A − 1 = 5.31 so λ₇ = trA − comp < 1, with G pivots
healthy (well-conditioned). Report the scan; pick best certified λ₇ < 1."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import arb_fejer_cert2 as ac
ctx.prec = 220
ac.W = arb('1.32')/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W; Wa = ac.W
trA_op = arb('6.3102004477213357703121798388447592044423984186896')

def sinc_sq(y):
    ay=abs(y); ymax=float(ay.mid())+float(ay.rad())
    if ymax<0.5:
        y2=y*y; r6=ymax**6/45+1e-30
        return acb(1)-y2/3+2*y2*y2/45+acb(arb(0,r6),arb(0,r6))
    return (y.sin()/y)**2
def Fev(c,r): return sinc_sq(Wa*(r-acb(c)))+sinc_sq(Wa*(r+acb(c)))

def compression(CEN):
    n=len(CEN)
    splits=sorted(set([14.0,44.0]+[float(c) for c in CEN if 14<c<44]))
    def AF(j,k):
        I=acb(0)
        for a,b in zip(splits[:-1],splits[1:]):
            I=I+acb.integral(lambda r,_: Fev(CEN[j],r)*Fev(CEN[k],r), acb(a),acb(b))
        return 2*I.real
    A=[[None]*n for _ in range(n)]
    for j in range(n):
        for k in range(j,n): A[j][k]=A[k][j]=AF(j,k)
    G=[[2*arb.pi()*ac.gpair_closed_arb(arb(CEN[j]),arb(CEN[k]),arb(0)) for k in range(n)] for j in range(n)]
    d=[None]*n; L=[[arb(0)]*n for _ in range(n)]
    for i in range(n):
        sd=G[i][i]
        for k in range(i): sd-=L[i][k]*L[i][k]*d[k]
        if not (sd>0): return None,None
        d[i]=sd
        for j2 in range(i+1,n):
            v=G[j2][i]
            for k in range(i): v-=L[j2][k]*L[i][k]*d[k]
            L[j2][i]=v/d[i]
    def solve(rhs):
        y=[None]*n
        for i in range(n):
            s=rhs[i]
            for k in range(i): s-=L[i][k]*y[k]
            y[i]=s
        for i in range(n): y[i]=y[i]/d[i]
        z=[None]*n
        for i in reversed(range(n)):
            s=y[i]
            for k in range(i+1,n): s-=L[k][i]*z[k]
            z[i]=s
        return z
    comp=arb(0)
    for col in range(n):
        x=solve([A[i][col] for i in range(n)]); comp+=x[col]
    return comp, min(float(d[i].mid()) for i in range(n))

for k in (6,8,10,12,14):
    CEN=[round(14+30*i/(k-1),3) for i in range(k)]
    comp, minpiv = compression(CEN)
    if comp is None:
        print(f"k={k}: G not PD (over capacity)"); continue
    lam7=trA_op-comp
    ub=float(lam7.upper())
    print(f"k={k:2d} sp={30/(k-1):.2f}: Σ={float(comp.mid()):.4f}±{float(comp.rad()):.0e} | Gpiv_min={minpiv:.2e} | λ₇≤{ub:.4f} {'✓<1 → ‖f‖²≤%.3f·T'%(1/(1-ub)) if ub<1 else '(weak)'}", flush=True)
