#!/usr/bin/env python3
"""§900: does the ~3% slack floor (§899) HOLD or THIN as frac→1 (critical sampling) at the low-height corner?
Push frac 0.6→0.97 at [14,44] 4W=4 & 5, tracking cond(G) (trust only where moderate), μ, π_w, ratio π_w/μ, margin.
If ratio stays ≥~1.02 as frac→1 (while cond(G) still moderate) → uniform floor plausible (proof has room). If ratio
→1 → razor asymptotically zero-slack. Honest: mark configs where cond(G) too high to trust."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh, pinvh
def vm(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=6000; LAM=vm(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)
def run(a,b,fourW,N,dr=0.003):
    W=fourW/4; RG=np.arange(-(b+70),(b+70),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    gw=np.linalg.eigvalsh(Gm); cond=gw[-1]/max(gw[0],1e-300)
    if gw[0]<=0: return None
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    gw2,gv=eigh(Gm); B=gv@np.diag(1/np.sqrt(np.maximum(gw2,1e-300)))@gv.T
    Msw=B@(arch-prime)@B; polew=B@pole@B
    ev,V=eigh(Msw); nneg=int((ev<0).sum())
    if nneg<1: return (N,cond,0,None,None,None,None)
    mu=-float(ev[0]); w=V[:,0]; Mplus=V[:,1:]@np.diag(ev[1:])@V[:,1:].T
    pev,pvec=eigh(polew); p=float(pev[-1]); v=pvec[:,-1]; pi_w=p*float(v@w)**2
    margin=float(eigh(Msw+polew,eigvals_only=True)[0])
    return (N,cond,nneg,mu,pi_w,pi_w/mu if mu>0 else None,margin)
print("="*104)
print("§900 frac→1 at low-height corner: does π_w/μ floor hold or thin? (trust only cond(G) < ~1e9)")
print("="*104)
for (a,b,fourW) in [(14,44,4.0),(14,44,5.0)]:
    Ncap=(b-a)*2*(fourW/4)/np.pi
    print(f"\n[{a},{b}] 4W={fourW} Ncap={Ncap:.1f}:")
    print(f"  {'N':>3} {'frac':>5} {'cond(G)':>9} {'#neg':>4} {'μ':>11} {'π_w':>11} {'π_w/μ':>7} {'margin':>11} {'trust':>6}")
    for frac in [0.55,0.65,0.75,0.85,0.90,0.94,0.97]:
        N=max(4,int(round(frac*Ncap))); r=run(a,b,fourW,N)
        if r is None: continue
        N,cond,nn,mu,pw,ratio,margin=r
        trust="yes" if cond<1e9 else "NOISE"
        rs=f"{ratio:.3f}" if ratio else "  -  "
        mus=f"{mu:.3e}" if mu else "   -   "; pws=f"{pw:.3e}" if pw else "   -   "; mgs=f"{margin:+.3e}" if margin is not None else "  -  "
        print(f"  {N:>3} {frac:>5.2f} {cond:>9.1e} {nn:>4} {mus:>11} {pws:>11} {rs:>7} {mgs:>11} {trust:>6}")
print("="*104)
print("ratio holding ≥~1.02 while cond(G)<1e9 as frac→1 → uniform floor plausible. ratio→1.00 → razor zero-slack.")
print("="*104)
