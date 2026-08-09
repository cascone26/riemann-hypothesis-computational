#!/usr/bin/env python3
"""§863: the Szego FLUCTUATION bound — a potential SUFFICIENT condition for positivity. min eigenvalue =
mean(S) - deviation. If (mean_S - min_eig)/std_S is BOUNDED by a constant C across windows, then
min_eig >= mean_S - C*std_S, so positivity holds when mean_S > C*std_S (mean beats C*spread). Both sides
concrete: mean_S = arch_mean - prime_mean; std_S = symbol std. Measure the ratio across windows/4W/N."""
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
def analyze(a,b,W,frac,dr=0.004):
    fourW=4*W; RGRID=np.arange(-(b+120),(b+120),dr); OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi)
    Ncap=(b-a)*2*W/np.pi; N=max(4,round(frac*Ncap)); CEN=np.linspace(a,b,N)
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
    mineig=float(eigh(pole+arch-prime,G,eigvals_only=True)[0])
    # symbol stats over window (weighted by nothing - uniform over [a,b])
    Psym=np.zeros_like(RGRID)
    for u,ww in zip(ln,lw): Psym+=ww*np.cos(u*RGRID)
    Psym*=(1/np.pi); S=OMEGA/(2*np.pi)-Psym; win=(RGRID>=a)&(RGRID<=b); Sw=S[win]
    meanS=float(Sw.mean()); stdS=float(Sw.std()); minS=float(Sw.min())
    return N, mineig, meanS, stdS, minS
print("="*104)
print("§863 Szego fluctuation: is (mean_S - min_eig)/std_S bounded? (would give sufficient condition mean_S > C*std_S)")
print("="*104)
print(f"{'window':>12} {'4W':>4} {'frac':>5} {'N':>3} {'min_eig':>11} {'mean_S':>8} {'std_S':>7} {'min_S':>8} {'(mean-mineig)/std':>17} {'mean/std':>9}")
ratios=[]
for (a,b,fourW) in [(14,44,2.0),(14,44,3.0),(50,80,3.0),(50,80,4.0),(120,150,4.0),(300,330,4.0)]:
    for frac in [0.5,0.75,0.95]:
        W=fourW/4; N,me,meanS,stdS,minS=analyze(a,b,W,frac)
        r=(meanS-me)/stdS; ratios.append(r)
        print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {frac:>5.2f} {N:>3} {me:>+11.4e} {meanS:>8.4f} {stdS:>7.4f} {minS:>+8.4f} {r:>17.3f} {meanS/stdS:>9.3f}")
print("="*104)
print(f"(mean_S - min_eig)/std_S range: [{min(ratios):.3f}, {max(ratios):.3f}]  -- if bounded by C, positivity when mean_S > C*std_S")
print("compare to mean/std column: if mean/std > max ratio in robust rows -> certified positivity structure")
print("="*104)
