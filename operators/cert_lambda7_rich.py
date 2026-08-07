#!/usr/bin/env python3
"""Certified λ₇ (§806 Handle): rich ambient basis (17 probes [12,44]/2, apparatus-consistent A & G),
project onto float top-6 concentration eigenvectors V (17×6), certify the 6×6 pencil trace
Σ = tr((VᵀGV)⁻¹(VᵀAV)) via interval LDLᵀ solve (VᵀGV is 6×6, well-separated → conditioned; NO
orthonormalization). Ky-Fan: Σ (any 6-dim subspace) ≤ Σ_top6(true) ⟹ λ₇ ≤ tr A − Σ. Split
integration at centers (no sinc 0/0). ALL module state set (TW/c_tri)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import arb_fejer_cert2 as ac
ctx.prec = 220
ac.W = arb('1.32')/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W; Wa = ac.W
CEN = [12+2*i for i in range(17) if 12+2*i <= 44]; n = len(CEN)
trA_op = arb('6.3102004477213357703121798388447592044423984186896')

def sinc_sq(y):
    ay=abs(y); ymax=float(ay.mid())+float(ay.rad())
    if ymax<0.5:
        y2=y*y; r6=ymax**6/45+1e-30
        return acb(1)-y2/3+2*y2*y2/45+acb(arb(0,r6),arb(0,r6))
    return (y.sin()/y)**2
def Fev(c,r): return sinc_sq(Wa*(r-acb(c)))+sinc_sq(Wa*(r+acb(c)))
splits=sorted(set([14,44]+[c for c in CEN if 14<c<44]))
def AF(j,k):
    I=acb(0)
    for a,b in zip(splits[:-1],splits[1:]):
        I=I+acb.integral(lambda r,_: Fev(CEN[j],r)*Fev(CEN[k],r), acb(a),acb(b))
    return 2*I.real
print(f"A ({n}x{n}, split integ)...", flush=True)
A=[[None]*n for _ in range(n)]
for j in range(n):
    for k in range(j,n): A[j][k]=A[k][j]=AF(j,k)
    print(f"  A row {j}", flush=True)
G=[[2*arb.pi()*ac.gpair_closed_arb(arb(CEN[j]),arb(CEN[k]),arb(0)) for k in range(n)] for j in range(n)]
# float top-6 eigvecs of (A,G)
Am=np.array([[float(A[j][k].mid()) for k in range(n)] for j in range(n)])
Gm=np.array([[float(G[j][k].mid()) for k in range(n)] for j in range(n)])
gw,gv=np.linalg.eigh(Gm); keep=gw>gw.max()*1e-11
Bf=gv[:,keep]/np.sqrt(gw[keep])
ew,ev=np.linalg.eigh(Bf.T@Am@Bf); order=np.argsort(ew)[::-1]
K=6; V=Bf@ev[:,order[:K]]
print(f"float top-6 concentration eigenvalues: {[f'{ew[order[i]]:.4f}' for i in range(6)]}", flush=True)
Va=[[arb(float(V[i][j])) for j in range(K)] for i in range(n)]
def VtMV(M):
    out=[[arb(0)]*K for _ in range(K)]
    for a in range(K):
        for b in range(a,K):
            s=arb(0)
            for i in range(n):
                for j2 in range(n): s+=Va[i][a]*M[i][j2]*Va[j2][b]
            out[a][b]=out[b][a]=s
    return out
MA=VtMV(A); MG=VtMV(G)
# LDLᵀ of MG (6x6), solve MG X = MA, trace
d=[None]*K; L=[[arb(0)]*K for _ in range(K)]
for i in range(K):
    sd=MG[i][i]
    for k in range(i): sd-=L[i][k]*L[i][k]*d[k]
    assert sd>0; d[i]=sd
    for j2 in range(i+1,K):
        v=MG[j2][i]
        for k in range(i): v-=L[j2][k]*L[i][k]*d[k]
        L[j2][i]=v/d[i]
print(f"  VᵀGV pivots min {min(float(d[i].mid()) for i in range(K)):.3e}", flush=True)
def solve(rhs):
    y=[None]*K
    for i in range(K):
        s=rhs[i]
        for k in range(i): s-=L[i][k]*y[k]
        y[i]=s
    for i in range(K): y[i]=y[i]/d[i]
    z=[None]*K
    for i in reversed(range(K)):
        s=y[i]
        for k in range(i+1,K): s-=L[k][i]*z[k]
        z[i]=s
    return z
comp=arb(0)
for col in range(K):
    x=solve([MA[i][col] for i in range(K)]); comp+=x[col]
lam7=trA_op-comp
print(f"\ncertified Σ₆ (compression on float top-6 subspace, ≤ Σ_top6) = {float(comp.mid()):.5f} ± {float(comp.rad()):.1e}")
print(f"===== certified λ₇ ≤ tr A − Σ₆ ≤ {float(lam7.upper()):.5f} =====")
if float(lam7.upper())<1:
    ff=1/(1-float(lam7.upper()))
    print(f"→ ‖f‖² ≤ {ff:.4f}·T   (k=6, CERTIFIED) — feeds assembly")
