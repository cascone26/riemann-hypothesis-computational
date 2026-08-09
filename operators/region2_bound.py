#!/usr/bin/env python3
"""§919: can Region 2 (shallow corners, 34<t<T0) be closed by a computable corrected bound? Crude §875 bound
margin≥floor−pnorm is negative there (fails). Test whether margin ≥ floor − pnorm_eff, where pnorm_eff = the prime
form's max eigenvalue RESTRICTED to the near-arch-floor subspace (the misalignment correction), closes Region 2. Also
report π_w/μ. Region 2 configs: shallow corners at moderate height."""
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
    if np.linalg.eigvalsh(Gm)[0]<=1e-13: return None
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    gw,gv=eigh(Gm); B=gv@np.diag(1/np.sqrt(gw))@gv.T
    APw=B@(arch+pole)@B; primew=B@prime@B
    floor=float(eigh(APw,eigvals_only=True)[0]); pnorm=float(eigh(primew,eigvals_only=True)[-1])
    margin=float(eigh(APw-primew,eigvals_only=True)[0])
    ev,V=eigh(B@(arch-prime)@B); mu=-ev[0] if ev[0]<0 else None
    # misalignment-corrected: project prime onto the low-eigenspace of arch+pole (the near-floor directions)
    aev,aV=eigh(APw); klow=aV[:,aev<floor+ (pnorm-floor)*0.5]  # directions where arch+pole is below the midpoint
    if klow.shape[1]>0:
        pnorm_eff=float(eigh(klow.T@primew@klow,eigvals_only=True)[-1])
    else: pnorm_eff=pnorm
    return N,floor,pnorm,margin,mu,pnorm_eff,floor-pnorm,floor-pnorm_eff
print("="*104)
print("§919 Region 2 corrected bound: margin vs (floor−pnorm crude) vs (floor−pnorm_eff misalignment-corrected)")
print("="*104)
print(f"{'window':>12} {'4W':>4} {'floor':>7} {'pnorm':>7} {'margin':>11} {'μ':>10} {'floor−pnorm':>11} {'floor−pnorm_eff':>15}")
for (a,b,fourW) in [(24,48,4.0),(30,60,5.0),(40,70,6.0),(45,75,7.0),(50,80,7.0),(60,90,8.0),(34,58,6.0),(40,64,7.0)]:
    r=run(a,b,fourW)
    if r is None: continue
    N,floor,pnorm,margin,mu,pneff,cr,corr=r
    mus=f"{mu:.2e}" if mu else " PSD "
    print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {floor:>7.4f} {pnorm:>7.4f} {margin:>+11.3e} {mus:>10} {cr:>+11.4f} {corr:>+15.4f}")
print("="*104)
print("crude floor−pnorm < 0 (fails, §875). If floor−pnorm_eff > 0 (corrected bound closes) in Region 2 → provable there.")
print("="*104)
