#!/usr/bin/env python3
"""§926 (autonomous, tactic #1): density-preserving decomposition. prime = SMOOTH (PNT main term
(1/π)∫_0^{4W} e^{u/2}cos(ru)du) + FLUCTUATION. Is arch+pole ⪰ smooth_prime COMFORTABLE (unlike the razor)?
If yes, RH-content = fluctuation. Compute margin(arch+pole−smooth) vs true margin across the corner."""
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
def build(a,b,fourW,frac=0.7,dr=0.004):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+60),(b+60),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    # true prime
    m=LOGN<fourW; ln=LOGN[m]; lw=LAMW[m]; cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N):
        prime[j,:]=2.0*(((P[j][None,:]*P)@cu.T)*dr/(2*np.pi)*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    # SMOOTH prime: symbol Psm(r)=(1/pi)∫_0^{4W} e^{u/2}cos(ru)du ; smooth form_jk = 2*(1/2pi)∫ P_j P_k * pi*Psm? 
    # match normalization to true: true symbol P(r)=(1/pi)Σ lw cos(ln r); smooth replaces Σ lw g(ln) -> ∫_0^{4W} e^{u/2} g(u) du
    uu=np.arange(0.0,fourW,0.002)
    smsym=np.array([ (np.exp(uu/2)*np.cos(r*uu)).sum()*0.002 for r in RG ])*(1.0/np.pi)
    smprime=np.zeros((N,N))
    for j in range(N):
        smprime[j,:]=2.0*((P[j][None,:]*P)@smsym)*dr/(2*np.pi)*np.pi  # 2*(1/2pi)*pi = the same 2/(2) factor; align below
    # NOTE normalization: true prime uses 2*(1/2pi)∫P_jP_k Σlw cos = (1/pi)... ; smprime uses same with Σlw->∫e^{u/2}
    smprime=0.5*(smprime+smprime.T)
    mtrue=float(eigh(arch+pole-prime,Gm,eigvals_only=True)[0])
    msm=float(eigh(arch+pole-smprime,Gm,eigvals_only=True)[0])
    return N,mtrue,msm,float(eigh(smprime,Gm,eigvals_only=True)[-1]),float(eigh(prime,Gm,eigvals_only=True)[-1])
print("="*92)
print("§926 density-preserving: margin(arch+pole − SMOOTH prime) vs true margin. is smooth comfortable?")
print("="*92)
print(f"{'window':>12} {'4W':>4} {'N':>3} {'true margin':>12} {'smooth margin':>13} {'pnorm(smooth)':>13} {'pnorm(true)':>11}")
for (a,b) in [(14,44),(24,54),(30,60),(45,75)]:
    for fourW in [3.0,4.0,6.0]:
        N,mt,ms,pns,pnt=build(a,b,fourW)
        print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {N:>3} {mt:>+12.3e} {ms:>+13.4f} {pns:>13.4f} {pnt:>11.4f}")
print("="*92)
print("smooth margin ≫ 0 (comfortable) while true margin razor → RH-content is the FLUCTUATION (smooth−true).")
print("smooth margin also razor/negative → even the PNT main term is razor (arch can't dominate smooth prime either).")
print("="*92)
