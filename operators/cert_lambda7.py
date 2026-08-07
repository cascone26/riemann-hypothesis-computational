#!/usr/bin/env python3
"""Certified λ₇ upper bound on the 6-probe capacity span (§806 support): A[j,k]=2∫_[14,44]F_jF_k
with integration domain SPLIT at every probe center (no rectangle straddles y=0 → no sinc 0/0).
Ky-Fan: λ₇ ≤ tr A − compression₆(G-orthonormal 6-probe span). tr A = 6.31024 (certified §804-fix)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import arb_fejer_cert2 as ac
ctx.prec = 200
ac.W = arb('1.32')/4; Wa = ac.W
CEN = [14,20,26,32,38,44]; n = len(CEN)
trA = arb('6.3102004477213357703121798388447592044423984186896')  # §804 certified even-op trace

def sinc_sq(y):  # acb, robust: series when |y| small (whole ball), closed form else
    ay = abs(y); ymax = float(ay.mid())+float(ay.rad())
    if ymax < 0.5:
        y2=y*y; r6=ymax**6/45+1e-30
        return acb(1)-y2/3+2*y2*y2/45 + acb(arb(0,r6),arb(0,r6))
    return (y.sin()/y)**2
def Fev(c, r): return sinc_sq(Wa*(r-acb(c))) + sinc_sq(Wa*(r+acb(c)))
# split points = window ends + all centers inside
splits = sorted(set([14,44]+[c for c in CEN if 14<c<44]))
def integ_AF(j,k):
    I=acb(0)
    for a,b in zip(splits[:-1], splits[1:]):
        I=I+acb.integral(lambda r,_: Fev(CEN[j],r)*Fev(CEN[k],r), acb(a), acb(b))
    return 2*I.real
print("A-matrix (split integration)...", flush=True)
A=[[None]*n for _ in range(n)]
for j in range(n):
    for k in range(j,n):
        A[j][k]=A[k][j]=integ_AF(j,k)
    print(f"  A row {j} done (diag rad {float(A[j][j].rad()):.1e})", flush=True)
# G (L2 overlap, closed form)
G=[[2*arb.pi()*ac.gpair_closed_arb(arb(CEN[j]),arb(CEN[k]),arb(0)) for k in range(n)] for j in range(n)]
# G-orthonormalize float, certify compression trace tr(U^T A U) with Neumann correction for U^TGU≈I
Gm=np.array([[float(G[j][k].mid()) for k in range(n)] for j in range(n)])
w,v=np.linalg.eigh(Gm); U=v@np.diag(w**-0.5)@v.T
Ua=[[arb(float(U[i][j])) for j in range(n)] for i in range(n)]
def UtMU(M):
    out=[[arb(0)]*n for _ in range(n)]
    for a in range(n):
        for b in range(a,n):
            s=arb(0)
            for i in range(n):
                for j2 in range(n): s+=Ua[i][a]*M[i][j2]*Ua[j2][b]
            out[a][b]=out[b][a]=s
    return out
M1=UtMU(A); M0=UtMU(G)
trM1=sum((M1[i][i] for i in range(n)), arb(0))
Efro=arb(0)
for i in range(n):
    for j2 in range(n):
        d=M0[i][j2]-(arb(1) if i==j2 else arb(0)); Efro+=d*d
Efro=Efro.sqrt()
M1fro=arb(0)
for i in range(n):
    for j2 in range(n): M1fro+=M1[i][j2]*M1[i][j2]
M1fro=M1fro.sqrt()
comp6 = trM1 - M1fro*Efro/(arb(1)-Efro)     # certified lower bound on sum of top-6 eigenvalues
lam7_ub = trA - comp6
print(f"\ncertified compression₆ (≥ S₆) = {float(comp6.mid()):.5f} ± {float(comp6.rad()):.1e}")
print(f"===== certified λ₇ ≤ tr A − S₆ ≤ {float(lam7_ub.upper()):.5f} =====")
print(f"→ ||f||² ≤ T/(1−λ₇) ≤ {1/(1-float(lam7_ub.upper())):.4f}·T  (k=6 assembly input, CERTIFIED)")
