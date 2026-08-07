#!/usr/bin/env python3
"""Clean certified λ₇ bound (§806 → Handle): NO orthonormalization. G is PD (L² Gram) → interval
LDLᵀ solve G·X = A gives X = G⁻¹A; tr(X) = Σλᵢ(A,G) = compression₆ (sum of the 6 concentration
eigenvalues on the capacity span). Ky-Fan: λ₇ ≤ tr A_op − compression₆. Split integration at
probe centers (no sinc 0/0). tr A_op = 6.31020 certified."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import arb_fejer_cert2 as ac
ctx.prec = 220
ac.W = arb('1.32')/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W; Wa = ac.W  # ALL module state (Red's §788 landmine: gpair reads TW/c_tri)
CEN = [14,20,26,32,38,44]; n = len(CEN)
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
print("A (windowed, split integ)...", flush=True)
A=[[None]*n for _ in range(n)]
for j in range(n):
    for k in range(j,n): A[j][k]=A[k][j]=AF(j,k)
    print(f"  A row {j}", flush=True)
G=[[2*arb.pi()*ac.gpair_closed_arb(arb(CEN[j]),arb(CEN[k]),arb(0)) for k in range(n)] for j in range(n)]
# diagnostic: A diag should be <= G diag (concentration <= full-line); ratio in [0,1]
for i in range(n):
    print(f"  A[{i}{i}]={float(A[i][i].mid()):.3f}  G[{i}{i}]={float(G[i][i].mid()):.3f}  ratio {float((A[i][i]/G[i][i]).mid()):.3f}")
# LDLᵀ of G (PD), then solve G X = A column-by-column, tr(X)
d=[None]*n; L=[[arb(0)]*n for _ in range(n)]
for i in range(n):
    sd=G[i][i]
    for k in range(i): sd-=L[i][k]*L[i][k]*d[k]
    assert sd>0, f"G not PD at pivot {i}"
    d[i]=sd
    for j2 in range(i+1,n):
        v=G[j2][i]
        for k in range(i): v-=L[j2][k]*L[i][k]*d[k]
        L[j2][i]=v/d[i]
print(f"  G pivots: min {min(float(d[i].mid()) for i in range(n)):.4e} (G well-conditioned if >0)")
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
trX=arb(0)
for col in range(n):
    x=solve([A[i][col] for i in range(n)]); trX+=x[col]
lam7=trA_op-trX
print(f"\ncompression₆ = tr(G⁻¹A) = {float(trX.mid()):.5f} ± {float(trX.rad()):.1e}  (= Σ of 6 concentration λ, must be ≤ {float(trA_op.mid()):.3f})")
print(f"===== certified λ₇ ≤ tr A − Σ₆ ≤ {float(lam7.upper()):.5f} =====")
if float(lam7.upper())<1:
    print(f"→ ‖f‖² ≤ {1/(1-float(lam7.upper())):.4f}·T  — k=6 assembly input CERTIFIED")
