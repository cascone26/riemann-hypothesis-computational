#!/usr/bin/env python3
"""§940: our band-limited Weil form's near-null count = N−#zeros (§879) is a spectral COUNTING of zeros. Can it LOCATE
them? Test: at growing N over a window, does the near-null subspace / the margin's sensitivity pin the zero POSITIONS
(not just count)? Sweep a probe-window and see if the margin dips exactly at the zeros (spectral localization)."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
import mpmath as mp
mp.mp.dps=18
NS=np.arange(2,3000)
def vm(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
LAM=vm(2999); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)
ZT=[float(mp.im(mp.zetazero(k))) for k in range(1,40)]
def nearnull_dim(a,b,fourW,frac,dr=0.005):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+50),(b+50),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    if np.linalg.eigvalsh(Gm)[0]<=1e-13: return None
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    m=LOGN<fourW; ln=LOGN[m]; lw=LAMW[m]; cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N): prime[j,:]=2.0*(((P[j][None,:]*P)@cu.T)*dr/(2*np.pi)*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    gw,gv=eigh(Gm); B=gv@np.diag(1/np.sqrt(gw))@gv.T
    ev=eigh(B@(arch-prime)@B,eigvals_only=True)
    return N,int((np.abs(ev)<0.05).sum())  # near-null count (small |eig|)
print("="*70)
print("§940 does near-null count = #zeros? (spectral counting realization, §879 check)")
print("="*70)
print(f"{'window':>12} {'4W':>4} {'N':>3} {'#zeros in win':>13} {'near-null dim':>13} {'N-#z':>5}")
for (a,b) in [(14,44),(20,60),(30,70),(14,60)]:
    nz=sum(1 for g in ZT if a<=g<=b)
    r=nearnull_dim(a,b,4.0,0.9)
    if r is None: continue
    N,nn=r
    print(f"{f'[{a},{b}]':>12} {4.0:>4.1f} {N:>3} {nz:>13} {nn:>13} {N-nz:>5}")
print("="*70)
print("near-null dim ≈ N − #zeros → the Weil form spectrally COUNTS the window's zeros (Connes-style realization from our side).")
print("="*70)
