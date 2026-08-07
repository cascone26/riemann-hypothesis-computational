#!/usr/bin/env python3
"""OVERNIGHT m_E(W) MAP + m_E·D structural-law test (§811). Dense scan of the certified pencil
margin m_E on the top-6 concentration eigenspace, across widths, probe count scaled to capacity
(resolved). Logs 4W, capacity=trA, nprobes, min pencil pivot, m_E, D(4W), m_E·D. Tests whether
m_E·D≈const (the fixed-margin-budget / Hodge-index-shape law §810b). Robust: try/except per W,
flush per point, flags pencil degradation. Runs for hours."""
import sys, os, math, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
ctx.prec = 200
t0=14.0
def Dcorner(fourW): return math.exp(fourW)-(t0/(2*math.pi))*(fourW-math.log(t0/(2*math.pi))+1)

def trA_even(fourW):
    band=arb(fourW)/2; K0=band/arb.pi()
    I2=acb.integral(lambda t,_: ((band*2*t).sin()/(arb.pi()*2*t)), acb(14), acb(44)).real
    return float((K0*30+I2).mid())

def mE_at(fourW, nprobes, K=6):
    ac.W=arb(fourW)/4; ac.TW=2*ac.W; ac.c_tri=arb.pi()/ac.W; Wa=ac.W
    CEN=[round(14+30*i/(nprobes-1),3) for i in range(nprobes)]; n=len(CEN); cs=[arb(c) for c in CEN]
    Q=[[None]*n for _ in range(n)]; G=[[None]*n for _ in range(n)]
    def sinc_sq(y):
        ay=abs(y); ymax=float(ay.mid())+float(ay.rad())
        if ymax<0.5:
            y2=y*y; return acb(1)-y2/3+2*y2*y2/45+acb(arb(0,ymax**6/45+1e-30),arb(0,ymax**6/45+1e-30))
        return (y.sin()/y)**2
    def Fev(c,r): return sinc_sq(Wa*(r-acb(c)))+sinc_sq(Wa*(r+acb(c)))
    splits=sorted(set([14.0,44.0]+[c for c in CEN if 14<c<44]))
    A=[[None]*n for _ in range(n)]
    for j in range(n):
        for k in range(j,n):
            Q[j][k]=Q[k][j]=ac.pole_entry(cs[j],cs[k])+zc.arch_entry_v3(cs[j],cs[k],arb(400))-ac.prime_entry(cs[j],cs[k])
            G[j][k]=G[k][j]=2*arb.pi()*ac.gpair_closed_arb(cs[j],cs[k],arb(0))
            I=acb(0)
            fine=[]
            for a,b in zip(splits[:-1],splits[1:]):
                m=max(2,int((b-a)))
                for tt in range(m): fine.append((a+(b-a)*tt/m,a+(b-a)*(tt+1)/m))
            for a,b in fine: I=I+acb.integral(lambda r,_: Fev(CEN[j],r)*Fev(CEN[k],r), acb(a),acb(b))
            A[j][k]=A[k][j]=2*I.real
    Am=np.array([[float(A[j][k].mid()) for k in range(n)] for j in range(n)])
    Gm=np.array([[float(G[j][k].mid()) for k in range(n)] for j in range(n)])
    gw,gv=np.linalg.eigh(Gm)
    if gw.min()<=0: return None,None  # G not PD
    Bf=gv/np.sqrt(gw); ew,ev=np.linalg.eigh(Bf.T@Am@Bf); order=np.argsort(ew)[::-1]
    V=Bf@ev[:,order[:K]]; Va=[[arb(float(V[i][j])) for j in range(K)] for i in range(n)]
    def VtMV(M):
        out=[[arb(0)]*K for _ in range(K)]
        for a in range(K):
            for b in range(a,K):
                s=arb(0)
                for i in range(n):
                    for j2 in range(n): s+=Va[i][a]*M[i][j2]*Va[j2][b]
                out[a][b]=out[b][a]=s
        return out
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
    if not ok(arb(0)): return None, float(gw.min())
    lo,hi=0.0,2.0
    for _ in range(42):
        mid=(lo+hi)/2
        if ok(arb(mid)): lo=mid
        else: hi=mid
    return lo, float(gw.min())

RESULTS="mE_overnight.txt"
open(RESULTS,"a").write("# 4W capacity nprobes Gpivmin m_E D_corner mE*D\n")
grid=[1.10,1.15,1.20,1.25,1.32,1.38,1.45,1.52,1.60,1.68,1.75,1.85,2.00,2.20,2.50,2.80,3.00]
for fourW in grid:
    try:
        cap=trA_even(fourW)
        npr=max(8, math.ceil(cap)+1)
        mE,gpiv=mE_at(fourW, npr)
        D=Dcorner(fourW)
        if mE is None:
            line=f"{fourW:.2f} {cap:.3f} {npr} {gpiv if gpiv else 0:.2e} DEGENERATE(pencil_or_G) - {D:.3f} -"
        else:
            line=f"{fourW:.2f} {cap:.3f} {npr} {gpiv:.2e} {mE:.6f} {D:.3f} {mE*D:.5f}"
        open(RESULTS,"a").write(line+"\n"); print(line, flush=True)
    except Exception as e:
        open(RESULTS,"a").write(f"{fourW:.2f} ERROR {type(e).__name__}\n"); print(f"{fourW:.2f} ERROR {e}", flush=True)
print("OVERNIGHT SCAN COMPLETE", flush=True)
