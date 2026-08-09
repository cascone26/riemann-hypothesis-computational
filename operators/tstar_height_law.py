#!/usr/bin/env python3
"""§857: does the prime-scaling safety margin t*(height) have a clean law in the robust regime?
t* grows with height (§856: 1.11,1.77,1.89). Sweep height at fixed 4W=4, N=13 (robust, N<#zeros for
height>=~60), compute t*, and see if it grows like the floor (~ln t), plateaus, or has another law."""
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
def mats(a,b,W,N,dr=0.004):
    fourW=4*W; RGRID=np.arange(-(b+120),(b+120),dr)
    OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi); CEN=list(np.linspace(a,b,N))
    def fej(x):
        y=W*x; out=np.ones_like(y); nz=np.abs(y)>1e-12; out[nz]=(np.sin(y[nz])/y[nz])**2; return out
    P=np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G=(P@P.T)*dr
    arch=(P*OMEGA[None,:])@P.T*dr/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cosu=np.cos(np.outer(ln,RGRID)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cosu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    return pole+arch,prime,G
def tstar(AP,prime,G):
    def ok(t): return float(eigh(AP-t*prime,G,eigvals_only=True)[0])>=0
    if not ok(0.0): return 0.0
    lo,hi=0.0,1.0
    while ok(hi) and hi<1e6: lo,hi=hi,hi*2
    for _ in range(46):
        mid=(lo+hi)/2
        if ok(mid): lo=mid
        else: hi=mid
    return lo
print("="*80); print("§857 t*(height) law, fixed 4W=4, N=13 (robust regime)"); print("="*80)
print(f"{'height':>8} {'floor':>8} {'t*':>9} {'t*/floor':>9} {'ln(t/2pi)':>10}")
ts=[]; hs=[]
for a in [60,120,250,500,1000,2000]:
    b=a+30; W=1.0; mid=(a+b)/2; AP,prime,G=mats(a,b,W,13)
    fl=float(eigh(AP,G,eigvals_only=True)[0]); t=tstar(AP,prime,G)
    ts.append(t); hs.append(mid)
    print(f"{mid:>8.0f} {fl:>8.4f} {t:>9.4f} {t/fl:>9.3f} {np.log(mid/(2*np.pi)):>10.3f}",flush=True)
ts=np.array(ts); hs=np.array(hs)
b1,b0=np.polyfit(np.log(hs),ts,1)
print("="*80); print(f"t* ~ {b1:.4f}*ln(height) + {b0:.4f}  (linear-in-ln fit)")
print(f"t* range {ts.min():.3f}-{ts.max():.3f}; {'PLATEAUS' if ts.max()-ts[-2]<0.05 else 'still growing'}")
print("="*80)
