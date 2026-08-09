#!/usr/bin/env python3
"""§914: SUB-LEMMA C — is the pole–worst-direction alignment (v·w)² bounded BELOW? v = pole direction (explicit,
= whitened f(i/2)), w = worst direction (M_σ neg eigvec). The corner criterion π_w=p(v·w)²≥μ; p explicit (§907), so
a LOWER bound on (v·w)² is the missing piece. Map (v·w)² across the bottom-razor region (h22-34, all 4W) and test:
is it bounded below? does it correlate with an explicit quantity (the pole overlap with the coincidence subspace)?"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
def vm(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=3200; LAM=vm(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)
def run(a,b,fourW,frac=0.85,dr=0.005):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+50),(b+50),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    gw=np.linalg.eigvalsh(Gm)
    if gw[0]<=1e-13: return None
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    gw2,gv=eigh(Gm); B=gv@np.diag(1/np.sqrt(gw2))@gv.T
    Msw=B@(arch-prime)@B; polew=B@pole@B; ev,V=eigh(Msw)
    if ev[0]>=0: return (N,None,None,None,None)
    mu=-ev[0]; w=V[:,0]; pev,pvec=eigh(polew); p=pev[-1]; v=pvec[:,-1]
    # overlap of pole direction v with the near-null (coincidence) subspace S (eigs |.|<1e-3 of M_σ)
    nn=np.abs(ev)<1e-3; US=V[:,nn]; v_in_S=np.linalg.norm(US.T@v)**2
    return (N,mu,(v@w)**2,v_in_S,p*(v@w)**2/mu)
print("="*92)
print("§914 sub-lemma C: is (v·w)² bounded below? (bottom-razor region, all 4W). also v's overlap w/ coincidence subspace")
print("="*92)
print(f"{'window':>12} {'4W':>4} {'N':>3} {'μ':>10} {'(v·w)²':>8} {'v in S':>7} {'π_w/μ':>7}")
allvw=[]
for fourW in [4.0,5.0,6.0,7.0,8.0]:
    for (a,b) in [(14,30),(14,44),(20,44),(24,48),(30,54)]:
        r=run(a,b,fourW)
        if r is None: continue
        N,mu,vw2,vinS,rt=r
        if mu is None: continue
        allvw.append(vw2)
        print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {N:>3} {mu:>10.3e} {vw2:>8.4f} {vinS:>7.3f} {rt:>7.3f}")
print("="*92)
if allvw:
    print(f"(v·w)² over {len(allvw)} bottom-razor corners: min={min(allvw):.4f} max={max(allvw):.4f} mean={np.mean(allvw):.4f}")
    print("If (v·w)² bounded below (min ≫ 0) AND 'v in S'≈(v·w)² → pole direction lives in coincidence subspace → boundable.")
print("="*92)
