"""Confirm the 15,25,35,45 basis positive margin is grid-converged (not a truncation artifact),
by widening/refining the r-grid and prime cutoff. If min-eig stabilizes at ~+0.048 => real,
certifiable margin. Also re-check 'wide' basis is a truncation artifact that heals with wider grid."""
import numpy as np, mpmath as mp
mp.mp.dps=25
from sympy import primerange
t=1.0
def phihat(a,r): return t*np.sqrt(2*np.pi)/2*(np.exp(-t**2*(r-a)**2/2)+np.exp(-t**2*(r+a)**2/2))
def phihat_c(a,r): return t*mp.sqrt(2*mp.pi)/2*(mp.e**(-t**2*(r-a)**2/2)+mp.e**(-t**2*(r+a)**2/2))
def g_cross_cf(aj,ak,y):
    s=aj+ak; d=aj-ak
    return np.exp(-y**2/4)*(np.sqrt(np.pi)/2)*(np.exp(-s**2/4)*np.cos(d*y/2)+np.exp(-d**2/4)*np.cos(s*y/2))
def primepowers(pmax,Lmax):
    pp=[]
    for p in primerange(2,pmax):
        lp=np.log(p); m=1
        while m*lp<Lmax: pp.append((m*lp, lp/np.sqrt(p**m))); m+=1
    return np.array([L for L,w in pp]),np.array([w for L,w in pp])
def buildQ(avals,rmax,rn,pmax,Lmax):
    Ls,Ws=primepowers(pmax,Lmax)
    rs=np.linspace(-rmax,rmax,rn)
    Om=np.array([float(mp.re(mp.digamma(mp.mpf('0.25')+1j*r/2))-mp.log(mp.pi)) for r in rs])
    nb=len(avals); Q=np.zeros((nb,nb))
    for j in range(nb):
        for k in range(nb):
            arch=np.trapezoid(phihat(avals[j],rs)*phihat(avals[k],rs)*Om,rs)/(2*np.pi)
            pole=2*float(mp.re(phihat_c(avals[j],mp.mpf('0.5')*1j)*phihat_c(avals[k],mp.mpf('0.5')*1j)))
            prime=2*np.sum(Ws*g_cross_cf(avals[j],avals[k],Ls))
            Q[j,k]=arch+pole-prime
    return (Q+Q.T)/2
av=[15.,25.,35.,45.]
print("basis 15,25,35,45 — grid convergence of min-eig:")
for rmax,rn,pmax,Lmax in [(60,9000,100000,25),(90,15000,100000,28),(120,24000,300000,30),(120,40000,300000,32)]:
    ev=np.linalg.eigvalsh(buildQ(av,rmax,rn,pmax,Lmax))
    print(f"  rmax={rmax} rn={rn} pmax={pmax} Lmax={Lmax}: min-eig={ev[0]:+.6e}  ev={np.array2string(ev,precision=4)}")
print("\n'wide' 20,40,60,80 healing with wider grid (should climb toward >=0):")
for rmax,rn in [(60,9000),(100,18000),(140,26000)]:
    ev=np.linalg.eigvalsh(buildQ([20.,40.,60.,80.],rmax,rn,300000,32))
    print(f"  rmax={rmax}: min-eig={ev[0]:+.6e}")
