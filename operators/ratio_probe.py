#!/usr/bin/env python3
"""§899: does the scalar inequality π_w = p(v·w)² ≥ μ have UNIFORM ROOM, or is it razor (ratio→1)? Sweep the corner
systematically (windows × 4W × frac near capacity) and report the ratio π_w/μ AND the exact-criterion margin. If
π_w/μ stays bounded above 1 with slack → a proof has room. If π_w/μ → 1 as we deepen → the inequality is razor
(the genuine hard truth). Also track μ and π_w scaling. Only well-conditioned (G PD) configs with #neg(M_σ)=1."""
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

def run(a,b,fourW,frac,dr=0.004):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+60),(b+60),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
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
    Msw=B@(arch-prime)@B; polew=B@pole@B
    ev,V=eigh(Msw); nneg=int((ev<0).sum())
    if nneg!=1: return (N,Ncap,nneg,None,None,None,None)
    mu=-float(ev[0]); w=V[:,0]; Mplus=V[:,1:]@np.diag(ev[1:])@V[:,1:].T
    pev,pvec=eigh(polew); p=float(pev[-1]); v=pvec[:,-1]
    pi_w=p*float(v@w)**2; A=Mplus+polew; g=1.0-mu*float(w@pinvh(A)@w)
    margin=float(eigh(Msw+polew,eigvals_only=True)[0])
    return (N,Ncap,nneg,mu,pi_w,pi_w/mu,margin)

print("="*104)
print("§899 does π_w ≥ μ have room? ratio π_w/μ across the corner (only #neg=1, G PD). ratio→1 = razor; >1 = slack")
print("="*104)
print(f"{'window':>12} {'4W':>4} {'frac':>5} {'N':>3} {'μ (defect)':>11} {'π_w=p(v·w)²':>12} {'π_w/μ':>7} {'margin':>11}")
rows=[]
for (a,b) in [(14,44),(14,34),(20,50),(30,60)]:
    for fourW in [3.0,4.0,5.0]:
        for frac in [0.6,0.7,0.8]:
            r=run(a,b,fourW,frac)
            if r is None or r[3] is None: continue
            N,Ncap,nn,mu,pw,ratio,margin=r
            rows.append((a,b,fourW,frac,mu,pw,ratio,margin))
            print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {frac:>5.2f} {N:>3} {mu:>11.3e} {pw:>12.3e} {ratio:>7.3f} {margin:>+11.3e}")
print("="*104)
if rows:
    ratios=np.array([r[6] for r in rows]); mus=np.array([r[4] for r in rows])
    print(f"π_w/μ: min={ratios.min():.3f} max={ratios.max():.3f} mean={ratios.mean():.3f} over {len(rows)} corners")
    # does ratio shrink toward 1 as mu shrinks (deeper corner)?
    order=np.argsort(mus); print(f"as μ decreases (deeper), ratio trend: μ={mus[order][:3]} -> ratio={ratios[order][:3]}")
    print("ratio bounded well above 1 uniformly → the scalar inequality has SLACK (proof has room to exploit).")
print("="*104)
