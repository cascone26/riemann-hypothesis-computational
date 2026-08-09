#!/usr/bin/env python3
"""§875: the rigorous sufficient condition floor(t) > pnorm(4W) => margin > 0. Verify: margin >= floor - pnorm
(operator bound), and floor-pnorm>0 covers the robust region (positivity rigorous), while the corner has
floor<pnorm (margin>0 only by misalignment, not covered). floor=lmin(arch+pole,G), pnorm=lmax(prime,G)."""
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
def parts(a,b,W,N,dr=0.004):
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
    floor=float(eigh(AP,G,eigvals_only=True)[0]); pnorm=float(eigh(prime,G,eigvals_only=True)[-1])
    margin=float(eigh(AP-prime,G,eigvals_only=True)[0])
    return floor,pnorm,margin
print("="*94)
print("§875 sufficient condition: floor(t) > pnorm(4W) => margin>0 (rigorous). margin >= floor-pnorm always.")
print("="*94)
print(f"{'window':>12} {'4W':>4} {'N':>3} {'floor':>8} {'pnorm':>8} {'floor-pnorm':>12} {'margin':>12} {'m>=f-p?':>8} {'suff?':>6}")
for (a,b,fourW,frac) in [(50,80,3.0,0.6),(120,150,4.0,0.6),(300,330,5.0,0.6),
                          (14,44,2.0,0.6),(14,44,3.0,0.7),(14,44,4.0,0.7),(60,90,4.0,0.7)]:
    W=fourW/4; L=b-a; N=max(4,round(frac*L*2*W/np.pi))
    fl,pn,mg=parts(a,b,W,N); fmp=fl-pn
    ok = mg >= fmp-1e-9
    suff = "YES pos" if fmp>0 else "corner"
    print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {N:>3} {fl:>8.4f} {pn:>8.4f} {fmp:>+12.4f} {mg:>+12.4e} {str(ok):>8} {suff:>6}")
print("="*94)
print("floor-pnorm>0 (robust region) => margin>0 RIGOROUSLY (unconditional sufficient condition).")
print("corner: floor<pnorm, margin>0 only via misalignment (not covered by this bound) = the hard part.")
