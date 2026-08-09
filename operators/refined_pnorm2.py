#!/usr/bin/env python3
"""§878: refined corner bound. Since v_p (prime top) is ~orthogonal to v_a (arch bottom, §835), the rank-1
downdate of the top prime eigenvalue doesn't lower the floor direction => margin ~>= floor - pnorm_2 (2nd prime
eigenvalue). Test: is floor > pnorm_2 at the CORNER (=> margin>0 by the refinement)? And does margin >= floor-pnorm_2?"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
def vonmangoldt(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=3000; LAM=vonmangoldt(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)
def probe(a,b,W,N,dr=0.004):
    fourW=4*W; RGRID=np.arange(-(b+120),(b+120),dr); OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi)
    CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; out=np.ones_like(y); nz=np.abs(y)>1e-12; out[nz]=(np.sin(y[nz])/y[nz])**2; return out
    P=np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G=(P@P.T)*dr; arch=(P*OMEGA[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cosu=np.cos(np.outer(ln,RGRID)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cosu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp)); AP=pole+arch
    gw,gv=eigh(G); B=gv@np.diag(1/np.sqrt(gw))@gv.T
    APw=B@AP@B; primew=B@prime@B
    floor=float(eigh(APw,eigvals_only=True)[0])
    pw=eigh(primew,eigvals_only=True)  # ascending
    pnorm=pw[-1]; pnorm2=pw[-2]; pnorm3=pw[-3]
    margin=float(eigh(APw-primew,eigvals_only=True)[0])
    # how many prime eigenvalues exceed the floor?
    n_above=int((pw>floor).sum())
    return floor,pnorm,pnorm2,pnorm3,margin,n_above
print("="*104)
print("§878 refined corner bound: is floor > pnorm_2 (=> margin>0 via v_p perp v_a rank-1 downdate)?")
print("="*104)
print(f"{'window':>10} {'4W':>4} {'N':>3} {'floor':>7} {'pnorm':>7} {'pnorm2':>7} {'pnorm3':>7} {'margin':>11} {'floor>pn2?':>10} {'#pn>floor':>9}")
for (a,b,fourW,N) in [(14,44,3.0,10),(14,44,4.0,13),(14,44,3.5,12),(50,80,4.0,12),(30,60,4.0,13),(14,44,2.5,9)]:
    W=fourW/4; fl,pn,pn2,pn3,mg,nab=probe(a,b,W,N)
    print(f"{f'[{a},{b}]':>10} {fourW:>4.1f} {N:>3} {fl:>7.4f} {pn:>7.4f} {pn2:>7.4f} {pn3:>7.4f} {mg:>+11.3e} {str(fl>pn2):>10} {nab:>9}")
print("="*104)
print("if floor>pnorm_2 at the corner => refined bound margin ~>= floor-pnorm_2 > 0 could PROVE corner positivity.")
print("#pn>floor = how many prime eigenvalues exceed floor (if just 1, the rank-1 refinement handles it).")
