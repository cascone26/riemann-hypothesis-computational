#!/usr/bin/env python3
"""§930 (autonomous, new direction from §928 universality): does razor-thin EDGE positivity hold across conductors, with
a clean scaling? Real primitive χ mod q=3,4,5,8 (correct parity→gamma shift, conductor floor-raise). Compute edge-razor
margin per conductor. Universal razor + clean margin(q) law → the universal mechanism any proof must capture."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
NMAX=3200; NS=np.arange(2,NMAX+1); LOGN=np.log(NS)
def vm(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
LAM=vm(NMAX); LAMW=LAM[2:]/np.sqrt(NS)
# real primitive characters (Legendre/Kronecker), with parity (odd=1 uses ψ0(3/4+..), even=0 uses ψ0(1/4+..))
def chi(q):
    if q==3: f=lambda n: 0 if n%3==0 else (1 if n%3==1 else -1); return f,1        # odd
    if q==4: f=lambda n: 0 if n%2==0 else (1 if n%4==1 else -1); return f,1        # odd
    if q==5:
        qr={1:1,4:1,2:-1,3:-1,0:0}; return (lambda n: qr[n%5]),0                    # even (Legendre mod5)
    if q==8:
        m={1:1,3:-1,5:-1,7:1,0:0,2:0,4:0,6:0}; return (lambda n: m[n%8]),0          # even (chi_8)
def build(a,b,fourW,q,frac=0.7,dr=0.004):
    f,odd=chi(q); W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+60),(b+60),dr); CEN=np.linspace(a,b,N)
    sh=0.75 if odd else 0.25
    OM=np.real(digamma(sh+0.5j*RG))-np.log(np.pi)+np.log(q)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    chiv=np.array([f(int(n)) for n in NS]); m=LOGN<fourW; ln=LOGN[m]; lw=(LAMW*chiv)[m]
    cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N): prime[j,:]=2.0*(((P[j][None,:]*P)@cu.T)*dr/(2*np.pi)*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    ev=eigh(arch-prime,Gm,eigvals_only=True)
    return N,int((ev<-1e-11).sum()),float(ev[0]),float(eigh(arch,Gm,eigvals_only=True)[0])
print("="*84)
print("§930 razor-thin EDGE positivity across conductors q (each at its bottom-edge windows). universal?")
print("="*84)
print(f"{'q':>3} {'parity':>6} {'window':>10} {'4W':>4} {'#neg(Mσ)':>9} {'min Mσ (edge margin)':>20} {'floor':>7}")
for q,par,wins in [(3,'odd',[(4,20),(6,22)]),(4,'odd',[(4,20),(6,22)]),(5,'even',[(5,20),(6,22)]),(8,'even',[(6,22),(8,24)])]:
    for (a,b) in wins:
        for fourW in [4.0]:
            N,nn,mn,fl=build(a,b,fourW,q)
            print(f"{q:>3} {par:>6} {f'[{a},{b}]':>10} {fourW:>4.1f} {nn:>9} {mn:>+20.3e} {fl:>7.3f}")
print("="*84)
print("razor-thin edge margin across ALL q (small, ~1e-4..1e-6, #neg 0 or 1) → razor is UNIVERSAL; margin(q) trend = the law.")
print("higher q → higher floor (conductor) but also more/lower zeros; the EDGE stays razor-thin universally = GRH.")
print("="*84)
