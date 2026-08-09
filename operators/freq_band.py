#!/usr/bin/env python3
"""§929 (workflow-suggested): is prime cancellation DECOMPOSABLE across frequency bands? Partition primes by log n into
bands; compute λmin(arch+pole − prime_band) per band vs the razor-thin global. Comfortable per-band + razor global →
decomposable (interface-limited, attackable). Razor/negative per-band too → global-irreducible (= RH, no lever)."""
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
def run(a,b,fourW,frac=0.7,dr=0.004):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+60),(b+60),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    def primeband(lo,hi):
        m=(LOGN>=lo)&(LOGN<hi)&(LOGN<fourW); ln=LOGN[m]; lw=LAMW[m]
        if len(ln)==0: return np.zeros((N,N))
        cu=np.cos(np.outer(ln,RG)); pr=np.zeros((N,N))
        for j in range(N): pr[j,:]=2.0*(((P[j][None,:]*P)@cu.T)*dr/(2*np.pi)*lw[None,:]).sum(axis=1)
        return 0.5*(pr+pr.T)
    bands=[(k,k+1) for k in range(int(np.ceil(fourW)))]
    bandmin=[float(eigh(arch+pole-primeband(lo,hi),Gm,eigvals_only=True)[0]) for lo,hi in bands]
    full=float(eigh(arch+pole-primeband(0,fourW),Gm,eigvals_only=True)[0])
    return N,bands,bandmin,full
print("="*90)
print("§929 frequency-band decomposition: λmin(arch+pole − prime_band) per band vs global (razor)")
print("="*90)
for (a,b,fourW) in [(14,44,4.0),(24,54,4.0),(30,60,4.0)]:
    N,bands,bandmin,full=run(a,b,fourW)
    print(f"\n[{a},{b}] 4W={fourW} N={N}: GLOBAL margin={full:+.3e}")
    for (lo,hi),bm in zip(bands,bandmin):
        pk=[int(NS[i]) for i in range(len(NS)) if lo<=LOGN[i]<hi and LOGN[i]<fourW][:6]
        print(f"    band log n∈[{lo},{hi}) (n={pk}): λmin(arch+pole−band) = {bm:+.4f}")
print("="*90)
print("Comfortable per-band + razor global → cancellation is at INTERFACES (decomposable, attackable).")
print("Per-band ALSO negative/razor → cancellation is GLOBAL-irreducible (= RH, no structural lever).")
print("="*90)
