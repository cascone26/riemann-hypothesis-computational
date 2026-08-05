"""Scan bases for the Weil positivity form, closed-form cross-correlation (fast+exact-friendly).
First VERIFY the closed form vs the numeric grid, then scan for a positive-margin basis."""
import numpy as np, mpmath as mp
mp.mp.dps=25
from sympy import primerange
t=1.0
def phix(a,x): return np.exp(-x**2/(2*t**2))*np.cos(a*x)
def phihat(a,r): return t*np.sqrt(2*np.pi)/2*(np.exp(-t**2*(r-a)**2/2)+np.exp(-t**2*(r+a)**2/2))
def phihat_c(a,r): return t*mp.sqrt(2*mp.pi)/2*(mp.e**(-t**2*(r-a)**2/2)+mp.e**(-t**2*(r+a)**2/2))
# closed form g_{jk}(y) (t=1)
def g_cross_cf(aj,ak,y):
    s=aj+ak; d=aj-ak
    return np.exp(-y**2/4)*(np.sqrt(np.pi)/2)*(np.exp(-s**2/4)*np.cos(d*y/2)+np.exp(-d**2/4)*np.cos(s*y/2))
# verify closed form vs grid for one (j,k)
xg=np.linspace(-30,30,12000)
def g_grid(aj,ak,y): return np.trapezoid(phix(aj,xg)*phix(ak,xg-y),xg)
aj,ak=4.,8.
err=max(abs(g_grid(aj,ak,yy)-g_cross_cf(aj,ak,yy)) for yy in [0.3,1.0,2.5,5.0])
print(f"closed-form check (a=4,8): max|grid-cf|={err:.2e}  {'OK' if err<1e-6 else 'MISMATCH'}")
pp=[]
for p in primerange(2,100000):
    lp=np.log(p); m=1
    while m*lp<25: pp.append((m*lp, lp/np.sqrt(p**m))); m+=1
Ls=np.array([L for L,w in pp]); Ws=np.array([w for L,w in pp])
rs=np.linspace(-60,60,9000)
Om=np.array([float(mp.re(mp.digamma(mp.mpf('0.25')+1j*r/2))-mp.log(mp.pi)) for r in rs])
def buildQ(avals):
    nb=len(avals); Q=np.zeros((nb,nb))
    for j in range(nb):
        for k in range(nb):
            arch=np.trapezoid(phihat(avals[j],rs)*phihat(avals[k],rs)*Om,rs)/(2*np.pi)
            pole=2*float(mp.re(phihat_c(avals[j],mp.mpf('0.5')*1j)*phihat_c(avals[k],mp.mpf('0.5')*1j)))
            prime=2*np.sum(Ws*g_cross_cf(avals[j],avals[k],Ls))
            Q[j,k]=arch+pole-prime
    return (Q+Q.T)/2
for name,av in [("saturated 0,4,8,12",[0.,4.,8.,12.]),
                ("on-zeros 14,18,22,26",[14.,18.,22.,26.]),
                ("on-zeros 15,25,35,45",[15.,25.,35.,45.]),
                ("wide 20,40,60,80",[20.,40.,60.,80.]),
                ("mix 0,15,30,45",[0.,15.,30.,45.])]:
    ev=np.linalg.eigvalsh(buildQ(av))
    print(f"{name:22s} min-eig={ev[0]:+.5e}  max={ev[-1]:+.4e}")
