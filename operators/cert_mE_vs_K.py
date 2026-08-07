#!/usr/bin/env python3
"""m_E(K) at FIXED W (§812): the clean 'margin per class' test §811 conflated. At one W, one resolved
basis, vary the concentration-eigenspace dimension K=2..nprobes; m_E(K)=λ_min(Q,G) on top-K. As K
grows (edge modes enter), m_E declines. Test m_E·K ≈ const (fixed-total-margin / Hodge-index law
measured PROPERLY with growing K). No W-confound, no nprobes steps. Run at 3 widths for W-stability."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
ctx.prec=200
def run(fourW, nprobes):
    ac.W=arb(fourW)/4; ac.TW=2*ac.W; ac.c_tri=arb.pi()/ac.W; Wa=ac.W
    CEN=[round(14+30*i/(nprobes-1),3) for i in range(nprobes)]; n=len(CEN); cs=[arb(c) for c in CEN]
    def sinc_sq(y):
        ay=abs(y); ymax=float(ay.mid())+float(ay.rad())
        if ymax<0.5:
            y2=y*y; return acb(1)-y2/3+2*y2*y2/45+acb(arb(0,ymax**6/45+1e-30),arb(0,ymax**6/45+1e-30))
        return (y.sin()/y)**2
    def Fev(c,r): return sinc_sq(Wa*(r-acb(c)))+sinc_sq(Wa*(r+acb(c)))
    splits=sorted(set([14.0,44.0]+[c for c in CEN if 14<c<44]))
    Q=[[None]*n for _ in range(n)]; G=[[None]*n for _ in range(n)]; A=[[None]*n for _ in range(n)]
    for j in range(n):
        for k in range(j,n):
            Q[j][k]=Q[k][j]=ac.pole_entry(cs[j],cs[k])+zc.arch_entry_v3(cs[j],cs[k],arb(400))-ac.prime_entry(cs[j],cs[k])
            G[j][k]=G[k][j]=2*arb.pi()*ac.gpair_closed_arb(cs[j],cs[k],arb(0))
            I=acb(0)
            for a,b in zip(splits[:-1],splits[1:]):
                m=max(2,int((b-a)))
                for tt in range(m): I=I+acb.integral(lambda r,_: Fev(CEN[j],r)*Fev(CEN[k],r), acb(a+(b-a)*tt/m),acb(a+(b-a)*(tt+1)/m))
            A[j][k]=A[k][j]=2*I.real
    Am=np.array([[float(A[j][k].mid()) for k in range(n)] for j in range(n)])
    Gm=np.array([[float(G[j][k].mid()) for k in range(n)] for j in range(n)])
    gw,gv=np.linalg.eigh(Gm); Bf=gv/np.sqrt(gw)
    ew,evc=np.linalg.eigh(Bf.T@Am@Bf); order=np.argsort(ew)[::-1]
    conc=[ew[order[i]] for i in range(n)]
    print(f"\n=== 4W={fourW}, nprobes={nprobes}, concentration eigs: {[f'{c:.3f}' for c in conc]} ===", flush=True)
    print(f"  K  concentration_λK   m_E(top-K)   m_E·K", flush=True)
    for K in range(2, n+1):
        V=Bf@evc[:,order[:K]]; Va=[[arb(float(V[i][j])) for j in range(K)] for i in range(n)]
        def VtMV(M):
            o=[[arb(0)]*K for _ in range(K)]
            for a in range(K):
                for b in range(a,K):
                    s=arb(0)
                    for i in range(n):
                        for j2 in range(n): s+=Va[i][a]*M[i][j2]*Va[j2][b]
                    o[a][b]=o[b][a]=s
            return o
        QE=VtMV(Q); GE=VtMV(G)
        def ok(lam):
            d=[None]*K; L=[[arb(0)]*K for _ in range(K)]
            for i in range(K):
                sd=QE[i][i]-lam*GE[i][i]
                for k in range(i): sd-=L[i][k]*L[i][k]*d[k]
                if not (sd>0): return False
                d[i]=sd
                for j2 in range(i+1,K):
                    v=QE[j2][i]-lam*GE[j2][i]
                    for k in range(i): v-=L[j2][k]*L[i][k]*d[k]
                    L[j2][i]=v/d[i]
            return True
        if not ok(arb(0)): print(f"  {K}  {conc[K-1]:.4f}   NOT PD", flush=True); continue
        lo,hi=0.0,3.0
        for _ in range(42):
            mid=(lo+hi)/2
            if ok(arb(mid)): lo=mid
            else: hi=mid
        print(f"  {K}  {conc[K-1]:.4f}   {lo:.6f}   {lo*K:.5f}", flush=True)
for fourW,npr in [(1.30,8),(1.46,8),(1.62,8)]:
    run(fourW,npr)
print("\nmE_vs_K COMPLETE", flush=True)
