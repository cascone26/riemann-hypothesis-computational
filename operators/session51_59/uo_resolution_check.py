# Is the lambda^2=60 break real or under-resolution? Re-test at larger N. If margin(Z)->0 and
# conical-max returns as N grows, it was discretization (a=2.05 large-a regime, LEDGER §45).
import numpy as np, scipy.linalg
from scipy.special import digamma as _digamma
TWO_PI=2*np.pi; LOGPI=float(np.log(np.pi))
def omega(r): return np.real(_digamma(0.25+0.5j*r))-LOGPI
def ehat(r,n,L):
    s=TWO_PI*n/L; d=s-r; out=np.empty_like(d,dtype=complex); sm=np.abs(d)<1e-12
    out[sm]=np.sqrt(L)+0j; dd=d[~sm]; out[~sm]=(1/np.sqrt(L))*np.exp(0.5j*dd*L)*2*np.sin(0.5*dd*L)/dd
    return out
def build_AB(idx,L,G):
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
        lp=float(np.log(p)); j=1
        while p**j<=LSQ+1e-9: c.append((float(p**j),lp)); j+=1
    return c
def buildP(idx,cmb,L):
    D=len(idx); P=np.zeros((D,D))
    for i,n in enumerate(idx):
        for j in range(i,D):
            m=idx[j]; P[i,j]=P[j,i]=sum(Lam*k**-0.5*qk(n,m,np.log(k),L) for k,Lam in cmb)
    return P
LSQ=60; L=float(2*np.log(np.sqrt(LSQ))); BP=[2,3,5,7,11,13,17,19,23,29,31,37]
rng=np.random.RandomState(1)
print(f"lambda^2=60 (a=logsqrt60={L/2:.2f}, large-a regime). resolution scan:")
for N in [10,16,22]:
    idx=list(range(-N,N+1)); G=max(20000, 1500*N)
    AB=build_AB(idx,L,G)
    def marg(pr): 
        tau=AB-buildP(idx,comb(pr,LSQ),L); tau=(tau+tau.T)/2; return scipy.linalg.eigvalsh(tau)[0]
    m0=marg(BP); c=np.ones(len(BP))/np.sqrt(len(BP)); bd=0; K=40; eps=0.05
    for _ in range(K):
        v=rng.randn(len(BP)); v=v-(v@c)*c; v/=np.linalg.norm(v)
        if marg([p*np.exp(eps*vi) for p,vi in zip(BP,v)])-m0<0 and marg([p*np.exp(-eps*vi) for p,vi in zip(BP,v)])-m0<0: bd+=1
    print(f"  N={N:2d} (D={2*N+1}): margin(Z)={m0:+.3e}   both-down={bd}/{K}")
print("margin(Z)->0 & both-down->K as N grows => lambda^2=60 break was under-resolution (universal optimality intact).")
