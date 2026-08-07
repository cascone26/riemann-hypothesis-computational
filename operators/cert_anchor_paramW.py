#!/usr/bin/env python3
"""SELF-CONSISTENT anchor theorem (§808): ONE conditioned 8-probe basis. Q (zero-free Weil), G (L²
pencil), A (concentration, split-integration) — all apparatus-consistent (module state set). E =
top-6 concentration eigenspace of (A,G). Certify on the SAME E: m_E = Q-pencil-λmin on E, λ₇ (7th
concentration eig), τ_E (E's weighted off-window mass). Assemble certified ε₀. Fixes the §806-λ
basis-consistency gap (m_E now on the true E, not the probe span)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
import os
FOURW = os.environ.get('FOURW','1.32')
NPROBE = int(os.environ.get('NPROBE','8'))
ctx.prec = 200
ac.W = arb(FOURW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W; Wa = ac.W
CEN = [round(14+30*i/(NPROBE-1),3) for i in range(NPROBE)]; n = len(CEN)
trA_op = arb('6.3102004477213357703121798388447592044423984186896')

def sinc_sq(y):
    ay=abs(y); ymax=float(ay.mid())+float(ay.rad())
    if ymax<0.5:
        y2=y*y; r6=ymax**6/45+1e-30
        return acb(1)-y2/3+2*y2*y2/45+acb(arb(0,r6),arb(0,r6))
    return (y.sin()/y)**2
def Fev(c,r): return sinc_sq(Wa*(r-acb(c)))+sinc_sq(Wa*(r+acb(c)))
splits=sorted(set([14.0,44.0]+[c for c in CEN if 14<c<44]))
cs=[arb(c) for c in CEN]

print("Q (zero-free), G (pencil), A (concentration) on 8-probe basis...", flush=True)
Q=[[None]*n for _ in range(n)]; G=[[None]*n for _ in range(n)]; A=[[None]*n for _ in range(n)]
for j in range(n):
    for k in range(j,n):
        Q[j][k]=Q[k][j]=ac.pole_entry(cs[j],cs[k])+zc.arch_entry_v3(cs[j],cs[k],arb(400))-ac.prime_entry(cs[j],cs[k])
        G[j][k]=G[k][j]=2*arb.pi()*ac.gpair_closed_arb(cs[j],cs[k],arb(0))
        I=acb(0)
        fine=[]
        for a,b in zip(splits[:-1],splits[1:]):
            m=max(2,int((b-a)))          # ~unit sub-segments, bounded count (no adaptive hang)
            for t in range(m): fine.append((a+(b-a)*t/m, a+(b-a)*(t+1)/m))
        for a,b in fine:
            I=I+acb.integral(lambda r,_: Fev(CEN[j],r)*Fev(CEN[k],r), acb(a),acb(b))
        A[j][k]=A[k][j]=2*I.real
    print(f"  row {j}", flush=True)

# E = top-6 concentration eigenspace of (A,G), float
Am=np.array([[float(A[j][k].mid()) for k in range(n)] for j in range(n)])
Gm=np.array([[float(G[j][k].mid()) for k in range(n)] for j in range(n)])
gw,gv=np.linalg.eigh(Gm); Bf=gv/np.sqrt(gw)
ew,ev=np.linalg.eigh(Bf.T@Am@Bf); order=np.argsort(ew)[::-1]
print(f"concentration eigenvalues: {[f'{ew[order[i]]:.4f}' for i in range(n)]}", flush=True)
K = max(4, NPROBE-2); V=Bf@ev[:,order[:K]]
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
QE=VtMV(Q); GE=VtMV(A if False else G); AE=VtMV(A)
# m_E = Q-pencil-λmin on E  (bisect)
def okpair(P,R,lam):
    d=[None]*K; L=[[arb(0)]*K for _ in range(K)]
    for i in range(K):
        sd=P[i][i]-lam*R[i][i]
        for k in range(i): sd-=L[i][k]*L[i][k]*d[k]
        if not (sd>0): return False
        d[i]=sd
        for j2 in range(i+1,K):
            v=P[j2][i]-lam*R[j2][i]
            for k in range(i): v-=L[j2][k]*L[i][k]*d[k]
            L[j2][i]=v/d[i]
    return True
assert okpair(QE,GE,arb(0)), "Q not PD on E"
lo,hi=0.0,2.0
for _ in range(45):
    mid=(lo+hi)/2
    if okpair(QE,GE,arb(mid)): lo=mid
    else: hi=mid
m_E=lo
# λ₇ via Σ₆ = tr(GE⁻¹ AE) ≤ Σ_top6 ⟹ λ₇ ≤ trA − Σ₆ (Ky-Fan; GE PD conditioned)
d=[None]*K; L=[[arb(0)]*K for _ in range(K)]
for i in range(K):
    sd=GE[i][i]
    for k in range(i): sd-=L[i][k]*L[i][k]*d[k]
    d[i]=sd
    for j2 in range(i+1,K):
        v=GE[j2][i]
        for k in range(i): v-=L[j2][k]*L[i][k]*d[k]
        L[j2][i]=v/d[i]
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
S6=arb(0)
for col in range(K):
    x=solve([AE[i][col] for i in range(K)]); S6+=x[col]
lam7=trA_op-S6
print(f"\ncertified m_E (Q-pencil-λmin on TRUE E) = {m_E:.6f}")
print(f"certified Σ₆ = {float(S6.mid()):.5f}, λ₇ ≤ {float(lam7.upper()):.5f}")
np.save(f'anchor_Evecs_{FOURW}.npy', V); np.save(f'anchor_CEN_{FOURW}.npy', np.array(CEN))
open('anchor_mW.txt','a').write(f"{FOURW} nprobes={NPROBE} K={K} m_E={m_E:.6f} lam_next={float(lam7.upper()):.6f}\n")
print("saved E-vectors for τ_E cert")
