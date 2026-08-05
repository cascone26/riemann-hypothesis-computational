"""ORGAN 2b — hierarchy: trace the Weil-form min-eigenvalue across growing/denser bases toward the
LOW-FREQUENCY (saturated) regime. RH => min-eig >= 0 always; the QUESTION is how it approaches 0
(the 'barely true' saturation) as the basis enriches. Closed-form cross-corr; wide grid (no trunc)."""
import numpy as np, mpmath as mp
mp.mp.dps=25
from sympy import primerange
t=1.0
def phihat(a,r): return np.sqrt(2*np.pi)/2*(np.exp(-(r-a)**2/2)+np.exp(-(r+a)**2/2))
def phihat_c(a,r): return mp.sqrt(2*mp.pi)/2*(mp.e**(-(r-a)**2/2)+mp.e**(-(r+a)**2/2))
def g_cf(aj,ak,y):
    s=aj+ak; d=aj-ak
    return np.exp(-y**2/4)*(np.sqrt(np.pi)/2)*(np.exp(-s**2/4)*np.cos(d*y/2)+np.exp(-d**2/4)*np.cos(s*y/2))
pp=[]
for p in primerange(2,300000):
    lp=np.log(p); m=1
    while m*lp<32: pp.append((m*lp, lp/np.sqrt(p**m))); m+=1
Ls=np.array([L for L,w in pp]); Ws=np.array([w for L,w in pp])
rmax=160; rs=np.linspace(-rmax,rmax,32000)
Om=np.array([float(mp.re(mp.digamma(mp.mpf('0.25')+1j*r/2))-mp.log(mp.pi)) for r in rs])
def buildQ(avals):
    nb=len(avals); Q=np.zeros((nb,nb))
    ph=[phihat(a,rs) for a in avals]
    for j in range(nb):
        for k in range(j,nb):
            arch=np.trapezoid(ph[j]*ph[k]*Om,rs)/(2*np.pi)
            pole=2*float(mp.re(phihat_c(avals[j],mp.mpf('0.5')*1j)*phihat_c(avals[k],mp.mpf('0.5')*1j)))
            prime=2*np.sum(Ws*g_cf(avals[j],avals[k],Ls))
            Q[j,k]=Q[k,j]=arch+pole-prime
    return Q
print("HIERARCHY A: dense low-frequency comb {0,Δ,...} covering [0,A], Δ=2  (into the saturation):")
for A in [4,8,12,16,20,24,30,40]:
    av=list(np.arange(0,A+0.1,2.0))
    ev=np.linalg.eigvalsh(buildQ(av))
    print(f"  [0..{A}] Δ2  dim={len(av):2d}  min-eig={ev[0]:+.5e}  2nd={ev[1]:+.4e}  max={ev[-1]:+.3e}")
print("\nHIERARCHY B: pure low-freq {0,1,2,3,4,5} (smoothest modes — where positivity is tightest):")
for A in [1,2,3,4,5,6,8]:
    av=list(np.arange(0,A+0.1,1.0))
    ev=np.linalg.eigvalsh(buildQ(av))
    print(f"  {{0..{A}}}  dim={len(av):2d}  min-eig={ev[0]:+.5e}")
