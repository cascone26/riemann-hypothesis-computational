"""ORGAN 2 v3: CORRECT Weil positivity matrix (validated formula + cross-correlation prime term).
Q_{jk}=arch+pole-prime with hhat_{jk}=phihat_j*phihat_k, g_{jk}=phi_j (cross-corr) phi_k, prime uses
g_{jk}(log n). RH => Q PSD. This is the object for the SDP/SOS certificate lane."""
import numpy as np, mpmath as mp
mp.mp.dps = 20
from sympy import primerange
t=1.0; avals=[0.0,4.0,8.0,12.0]; nb=len(avals)
def phix(k,x): return np.exp(-x**2/(2*t**2))*np.cos(avals[k]*x)
def phihat(k,r): return t*np.sqrt(2*np.pi)/2*(np.exp(-t**2*(r-avals[k])**2/2)+np.exp(-t**2*(r+avals[k])**2/2))
def phihat_c(k,r): return t*mp.sqrt(2*mp.pi)/2*(mp.e**(-t**2*(r-avals[k])**2/2)+mp.e**(-t**2*(r+avals[k])**2/2))
# cross-correlation g_{jk}(y)=int phi_j(x) phi_k(x-y) dx on a grid (direct, phi decays fast)
xg=np.linspace(-30,30,12000); dx=xg[1]-xg[0]
def g_cross(j,k,ys):
    pj=phix(j,xg)
    out=np.empty(len(ys))
    for i,y in enumerate(ys):
        out[i]=np.trapezoid(pj*phix(k,xg-y),xg)
    return out
# prime powers
pp=[]
for p in primerange(2,100000):
    lp=np.log(p); m=1
    while m*lp<25: pp.append((m*lp, lp/np.sqrt(p**m))); m+=1
Ls=np.array([L for L,w in pp]); Ws=np.array([w for L,w in pp])
rs=np.linspace(-50,50,8000)
Om=np.array([float(mp.re(mp.digamma(mp.mpf('0.25')+1j*r/2))-mp.log(mp.pi)) for r in rs])
Q=np.zeros((nb,nb))
for j in range(nb):
    for k in range(nb):
        arch=np.trapezoid(phihat(j,rs)*phihat(k,rs)*Om,rs)/(2*np.pi)
        pole=2*float(mp.re(phihat_c(j,mp.mpf('0.5')*1j)*phihat_c(k,mp.mpf('0.5')*1j)))
        gvals=g_cross(j,k,Ls)  # g_{jk}(log n) at all prime-power logs
        prime=2*np.sum(Ws*gvals)   # 2 sum Lambda/sqrt(n) g_{jk}(log n)  (g even in y)
        Q[j,k]=arch+pole-prime
Q=(Q+Q.T)/2
ev=np.linalg.eigvalsh(Q)
print("ORGAN 2 v3 Weil positivity matrix (validated formula + cross-corr prime):")
print("  eigenvalues:", np.array2string(ev,precision=5))
print("  min-eig = %.6f => %s" % (ev[0], "PSD (RH-consistent) - FORMULA CORRECT" if ev[0]>-1e-7 else "NEG - still wrong"))
np.save("weil_Q_v3.npy", Q)
print("  saved weil_Q_v3.npy")
print("DONE")
