# §64: UNIVERSAL-OPTIMALITY test of the LP/energy lead. Is Z the conical margin-MAX across MANY
# bandwidths lambda^2 (the E8-signature: optimal for a whole family, not one functional)? If yes =>
# an interpolation/magic-function certificate is plausible; if it holds only at lambda^2=14 => not universal, lead dies.
import numpy as np, scipy.linalg
from scipy.special import digamma as _digamma
TWO_PI=2*np.pi; LOGPI=float(np.log(np.pi)); G=20000
def omega(r): return np.real(_digamma(0.25+0.5j*r))-LOGPI
def ehat(r,n,L):
    s=TWO_PI*n/L; d=s-r; out=np.empty_like(d,dtype=complex); sm=np.abs(d)<1e-12
    out[sm]=np.sqrt(L)+0j; dd=d[~sm]; out[~sm]=(1/np.sqrt(L))*np.exp(0.5j*dd*L)*2*np.sin(0.5*dd*L)/dd
    return out
def build_AB(idx,L):
    R=TWO_PI*max(abs(n) for n in idx)/L+80.0; grid=np.linspace(-R,R,G); h=grid[1]-grid[0]
    w=np.full(G,h); w[0]=w[-1]=0.5*h; Om=omega(grid); E=np.array([ehat(grid,n,L) for n in idx])
    A=((E*(w*Om))@E.conj().T).real/(2*np.pi)
    D=len(idx); Bm=np.zeros((D,D)); sh=np.sinh(L/4)**2
    for i,n in enumerate(idx):
        for j,m in enumerate(idx):
            num=L**2-16*np.pi**2*m*n; den=(L**2+16*np.pi**2*m**2)*(L**2+16*np.pi**2*n**2)
            Bm[i,j]=32*L*sh*num/den
    return A+Bm
def qk(n,m,y,L):
    if n==m: return 2*(1-abs(y)/L)*np.cos(TWO_PI*n*y/L)
    return (np.sin(TWO_PI*m*y/L)-np.sin(TWO_PI*n*y/L))/(np.pi*(n-m))
def comb(primes,LSQ):
    c=[]
    for p in primes:
        if p<=1: continue
        lp=float(np.log(p)); j=1
        while p**j<=LSQ+1e-9: c.append((float(p**j),lp)); j+=1
    return c
def buildP(idx,cmb,L):
    D=len(idx); P=np.zeros((D,D))
    for i,n in enumerate(idx):
        for j in range(i,D):
            m=idx[j]; P[i,j]=P[j,i]=sum(Lam*k**-0.5*qk(n,m,np.log(k),L) for k,Lam in cmb)
    return P
def primes_upto(x): return [p for p in [2,3,5,7,11,13,17,19,23,29,31,37] if p<=x]
rng=np.random.RandomState(1); N=10; idx=list(range(-N,N+1))
print(f"{'lam^2':>6} {'#primes':>7} {'margin(Z)':>13} {'both-down/60 (DENSITY-PRESERVING)':>34}")
for LSQ in [6,14,30,60]:
    L=float(2*np.log(np.sqrt(LSQ))); AB=build_AB(idx,L)
    BP=primes_upto(LSQ)
    def marg(pr): 
        tau=AB-buildP(idx,comb(pr,LSQ),L); tau=(tau+tau.T)/2; return scipy.linalg.eigvalsh(tau)[0]
    m0=marg(BP); c=np.ones(len(BP))/np.sqrt(len(BP))
    bd=0; K=60; eps=0.05
    for _ in range(K):
        v=rng.randn(len(BP)); v=v-(v@c)*c; nn=np.linalg.norm(v)
        if nn<1e-9: continue
        v/=nn
        if marg([p*np.exp(eps*vi) for p,vi in zip(BP,v)])-m0<0 and marg([p*np.exp(-eps*vi) for p,vi in zip(BP,v)])-m0<0: bd+=1
    print(f"{LSQ:>6} {len(BP):>7} {m0:>13.2e} {bd:>25}/60")
print("\nAll rows ~60/60 => Z is UNIVERSALLY the conical margin-max across bandwidths => interpolation certificate plausible.")
