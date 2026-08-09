#!/usr/bin/env python3
"""§887: does the near-null Weil eigenfunction φ* CONCENTRATE on the symbol's negative set {σ<0}? If the
margin-achieving band-limited φ* piles its |φ*|² mass on σ<0 (vs a generic band-limited fn), the corner is a
concentration/uncertainty problem: margin = floor − (negative mass φ* can capture). Compute φ* (min gen-eigvec
of Weil form vs G in Fejér basis), reconstruct φ*(r), measure P_neg(φ*)=∫_{σ<0}|φ*|²/∫|φ*|² vs P_neg(generic)."""
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

def run(a,b,fourW,N,dr=0.002):
    W=fourW/4; RGRID=np.arange(-(b+120),(b+120),dr)
    OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi); CEN=np.linspace(a,b,N)
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
    Q=arch+pole-prime
    w,V=eigh(Q,G); coef=V[:,0]  # min gen-eigenvector = φ*
    margin=w[0]
    # reconstruct φ*(r) on the window grid (use r>0 side, |φ*|^2)
    Pw=(1/np.pi)*(lw[:,None]*np.cos(np.outer(ln,RGRID))).sum(axis=0)
    sig=OMEGA/(2*np.pi)-Pw
    win=(RGRID>=a-3)&(RGRID<=b+3)
    phi=(coef[:,None]*P).sum(axis=0)  # φ*(r)
    dens=phi**2
    negmask=sig<0
    # concentration on {σ<0}, restricted to the window region where the mode lives
    reg=win
    tot=(dens[reg]).sum(); Pneg_star=(dens[reg & negmask]).sum()/tot
    # generic band-limited fn: average P_neg over the N basis functions
    Pneg_gen=np.mean([ (P[j][reg&negmask]**2).sum()/max((P[j][reg]**2).sum(),1e-30) for j in range(N)])
    # also: value ∫|φ*|² σ  (should be ~ margin*normalization) and where negative mass sits
    num=(dens[reg]*sig[reg]).sum()*dr; den=(dens[reg]).sum()*dr
    return margin,Pneg_star,Pneg_gen,num/den,sig[reg].min()

print("="*100)
print("§887 does φ* concentrate on {σ<0}?  P_neg=mass fraction on negative-symbol set  (φ* vs generic band fn)")
print("="*100)
print(f"{'window':>12} {'4W':>4} {'N':>3} {'margin':>11} {'P_neg(φ*)':>10} {'P_neg(gen)':>11} {'ratio':>6} {'<σ>_φ*':>9}")
for (a,b,fourW,N) in [(14,44,4.0,17),(14,44,3.0,13),(14,44,5.0,20),(60,90,4.0,30),(120,150,4.0,38)]:
    m,ps,pg,avg,smin=run(a,b,fourW,N)
    print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {N:>3} {m:>+11.3e} {ps:>10.3f} {pg:>11.3f} {ps/pg:>6.2f} {avg:>+9.4f}")
print("="*100)
print("P_neg(φ*) ≫ P_neg(generic) → φ* deliberately concentrates on σ<0 → corner IS a concentration/uncertainty")
print("bound: floor(positive mean) vs the negative mass a bandwidth-4W fn can capture. <σ>_φ* ≈ margin/‖φ*‖² > 0.")
print("="*100)
