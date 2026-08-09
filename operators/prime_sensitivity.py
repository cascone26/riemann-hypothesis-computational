#!/usr/bin/env python3
"""§860: prime-sensitivity — which primes are 'load-bearing' for Weil positivity? At a robust window
(margin ~0.1), perturb each prime's weight Lambda(n)/sqrt(n) by a factor and find the critical factor that
drives the margin to 0. Small critical perturbation = load-bearing prime. Maps the arithmetic structure of
positivity. Valid (real primes, full Weil form; only a single weight is scaled)."""
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
NMAX=3000; LAM=vonmangoldt(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW0=LAM[2:]/np.sqrt(NS)
a,b,W=14.0,44.0,0.5   # 4W=2.0, robust (margin ~0.1)
fourW=4*W; N=7; dr=0.004
RGRID=np.arange(-(b+120),(b+120),dr); OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi)
CEN=np.linspace(a,b,N)
def fej(x):
    y=W*x; out=np.ones_like(y); nz=np.abs(y)>1e-12; out[nz]=(np.sin(y[nz])/y[nz])**2; return out
P=np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G=(P@P.T)*dr
arch=(P*OMEGA[None,:])@P.T*dr/(2*np.pi)
def fc(c):
    y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
mask=LOGN<fourW; ln=LOGN[mask]; idxs=np.where(mask)[0]
cosu=np.cos(np.outer(ln,RGRID))
# precompute per-n prime matrices M_n so prime = sum_n w_n M_n
Mn=[]
for jj in range(len(ln)):
    col=cosu[jj]
    Mk=np.zeros((N,N))
    for j in range(N):
        Mk[j,:]=2.0*(((P[j][None,:]*P)@col)*dr/(2*np.pi))
    Mn.append(0.5*(Mk+Mk.T))
Mn=np.array(Mn)
w0=LAMW0[mask]
def margin(weights):
    prime=np.tensordot(weights,Mn,axes=1)
    return float(eigh(pole+arch-prime,G,eigvals_only=True)[0])
base=margin(w0)
print(f"§860 prime-sensitivity at [14,44] 4W=2.0 (base margin={base:+.4e})")
print(f"{'n':>4} {'weight':>9} {'crit_factor(margin->0)':>22} {'interp':>18}")
# for each prime-power n in range, scale its weight by f, find f where margin=0
nn=NS[mask]
for i in range(len(ln)):
    def m_of_f(f):
        w=w0.copy(); w[i]=w0[i]*f; return margin(w)
    # margin decreases as we ADD prime weight (prime is subtracted); find f>1 where margin->0, or f<0
    lo,hi=1.0,1.0; 
    # search upward: increasing this weight
    f=1.0; step=1.0; last=base; crit=None
    for f in np.arange(1.0,60.0,0.5):
        if m_of_f(f)<0: crit=f; break
    tag = f"LOAD-BEARING" if (crit and crit<5) else ("sensitive" if crit and crit<20 else "robust")
    cs = f"{crit:.2f}" if crit else ">60"
    print(f"{nn[i]:>4} {w0[i]:>9.4f} {cs:>22} {tag:>18}")
print("small crit_factor = perturbing THIS prime's weight easily breaks positivity (load-bearing).")
