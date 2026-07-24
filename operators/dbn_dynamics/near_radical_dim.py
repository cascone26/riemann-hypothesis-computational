# Connes Step 2, lemma L1 (LEDGER §100): what is the DIMENSION of the QW_lambda near-radical on [lambda^-1,lambda],
# and which prolate modes span it? Probe: for each even prolate mode h_{n,lambda}, compute the Weil Rayleigh quotient
# of E(h_{n,lambda}) and classify near-null (<~floor) vs bulk (O(1)). The float64 floor (~1e-7) does NOT defeat a
# COUNT (how many modes are small). If only a FEW low even modes are near-null, the near-radical is low-dim (L1),
# and the h_0,h_4 selection is by the xi-matching (k^=Xi), not by concentration alone.
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto
from numpy.polynomial.legendre import legval

def y2_matrix(Nleg):
    Y=np.zeros((Nleg+2,Nleg+2))
    for n in range(Nleg+1):
        if n+1<=Nleg+1: Y[n+1,n]=(n+1)/(2*n+1)*np.sqrt((2*n+1)/(2*n+3))
        if n-1>=0:      Y[n-1,n]=n/(2*n+1)*np.sqrt((2*n+1)/(2*n-1))
    return (Y[:Nleg,:Nleg])@(Y[:Nleg,:Nleg])
def prolate(lmbda,Nleg):
    c=2*np.pi*lmbda**2; n=np.arange(Nleg)
    A=np.diag((n*(n+1)).astype(float))+c**2*y2_matrix(Nleg); A=(A+A.T)/2
    w,V=np.linalg.eigh(A)
    par=np.array([(-1)**int(np.argmax(np.abs(V[:,k]))) for k in range(V.shape[1])])
    return V,np.where(par>0)[0],Nleg

def weil_rayleigh(Efunc_on_s, a, Lam, Nmax, M=6000):
    s=np.linspace(-a,a,M+1); h=s[1]-s[0]
    diffs=np.arange(-M,M+1)*h
    gvals=np.array([screw_g(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1); G=gvals[(idx[:,None]-idx[None,:])+M]
    k=Efunc_on_s(s); kp=np.gradient(k,s)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wk=kp*w; Q=float(Wk@G@Wk); n2=float(np.sum(k*k*w))
    return Q/n2 if n2>0 else np.nan

print("L1 probe: Weil Rayleigh of E(h_{n,lambda}) per even prolate mode n. near-null vs bulk => near-radical dimension.")
for lmbda in [1.6, 2.0]:
    a=np.log(lmbda); Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    Nleg=max(200,int(6*2*np.pi*lmbda**2)); V,ei,Nl=prolate(lmbda,Nleg)
    print(f"\n lambda={lmbda} (a={a:.3f}, ~4l^2={4*lmbda**2:.0f} concentrated modes):")
    print(f"   {'evenmode':>8} {'QW(E(h_n))/||.||^2':>18} {'class':>10}")
    for j,col in enumerate(ei[:7]):
        coeff=np.array([V[m,col]*np.sqrt((2*m+1)/2) for m in range(Nl)])
        def Ef(s, coeff=coeff):
            u=np.exp(s); nmax=int(np.floor(lmbda**2))+1; ssum=np.zeros_like(u)
            for nn in range(1,nmax):
                y=nn*u/lmbda; msk=np.abs(y)<=1; t=np.zeros_like(u)
                t[msk]=legval(y[msk],coeff); ssum+=t
            return np.sqrt(u)*ssum
        R=weil_rayleigh(Ef,a,Lam,Nmax)
        cls='near-null' if abs(R)<1e-5 else ('small' if abs(R)<1e-2 else 'BULK')
        print(f"   h_{{{2*j},l}}  {R:>18.4e} {cls:>10}",flush=True)
print("\nREAD: how many even modes are near-null? That = near-radical dim (L1). If it's small & low (h_0,h_2,h_4..),")
print("the radical is finite/low-dim and the h_0,h_4 SELECTION for k_lambda is the xi-matching (k^=Xi), not concentration.")
