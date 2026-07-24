# Connes Step 2 lemma (LEDGER §101), analytic CORE: is the QW_lambda-minimizer over a span of prolate modes,
# SUBJECT TO f(0)=fhat(0)=0, the k_lambda (h_0,h_4) direction? Mechanism (Connes/CMbook): E(f) in radical of GLOBAL
# Weil form (QW_inf(E f)=0 by Poisson), so QW_lambda(E f) = boundary terms of the truncation to [l^-1,l]; the two
# conditions f(0)=fhat(0)=0 kill the O(1) boundary terms, leaving O(1-chi). If the constrained minimizer over
# {h_0..h_10} collapses onto h_0,h_4, that IS the selection mechanism (and it's WELL-CONDITIONED: 6x6 in the prolate
# basis, immune to the float64 floor that kills the raw 1e-14 eigenvalue).
import numpy as np
from numpy.polynomial.legendre import legval
import sys; sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def y2_matrix(N):
    Y=np.zeros((N+2,N+2))
    for n in range(N+1):
        if n+1<=N+1: Y[n+1,n]=(n+1)/(2*n+1)*np.sqrt((2*n+1)/(2*n+3))
        if n-1>=0:   Y[n-1,n]=n/(2*n+1)*np.sqrt((2*n+1)/(2*n-1))
    return (Y[:N,:N])@(Y[:N,:N])
def prolate(lmbda,N):
    c=2*np.pi*lmbda**2; n=np.arange(N)
    A=np.diag((n*(n+1)).astype(float))+c**2*y2_matrix(N); A=(A+A.T)/2
    w,V=np.linalg.eigh(A); par=np.array([(-1)**int(np.argmax(np.abs(V[:,k]))) for k in range(V.shape[1])])
    return V,np.where(par>0)[0]

def analyze(lmbda, nmodes=6, M=6000):
    a=np.log(lmbda); Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000)); N=max(200,int(6*2*np.pi*lmbda**2))
    V,ei=prolate(lmbda,N); cols=ei[:nmodes]
    s=np.linspace(-a,a,M+1); h=s[1]-s[0]
    diffs=np.arange(-M,M+1)*h
    gv=np.array([screw_g(d,Lam,Nmax) for d in diffs]); idx=np.arange(M+1)
    Gk=gv[(idx[:,None]-idx[None,:])+M]
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    # E(h_{n,l}) on s-grid, and its derivative
    def Emap(col):
        coeff=np.array([V[m,col]*np.sqrt((2*m+1)/2) for m in range(N)])
        u=np.exp(s); nmax=int(np.floor(lmbda**2))+1; ssum=np.zeros_like(u)
        for nn in range(1,nmax):
            y=nn*u/lmbda; msk=np.abs(y)<=1; t=np.zeros_like(u)
            t[msk]=legval(y[msk],coeff); ssum+=t
        return np.sqrt(u)*ssum
    E=[Emap(c) for c in cols]; Ep=[np.gradient(e,s) for e in E]
    n=len(cols)
    Gw=np.zeros((n,n)); Mm=np.zeros((n,n))          # Weil Gram, mass Gram
    for i in range(n):
        for j in range(n):
            Gw[i,j]=(Ep[i]*w)@Gk@(Ep[j]*w)
            Mm[i,j]=np.sum(E[i]*E[j]*w)
    Gw=(Gw+Gw.T)/2; Mm=(Mm+Mm.T)/2
    # constraints on prolate coeffs c: f=sum c_i h_{col_i}. f(0)=sum c_i h_{col_i}(0); fhat(0)=int f = sum c_i * (sqrt2 * l * V[0,col_i])
    h0=np.array([legval(0.0,np.array([V[m,c]*np.sqrt((2*m+1)/2) for m in range(N)])) for c in cols])  # h_{n,l}(0)
    intg=np.array([np.sqrt(2)*lmbda*V[0,c] for c in cols])                                            # int h_{n,l}
    return Gw,Mm,h0,intg,cols

def constrained_min(Gw,Mm,C):
    """min c^T Gw c / c^T Mm c  s.t. C c = 0. Project onto ker(C)."""
    # orthonormal basis of ker(C)
    u,sv,vt=np.linalg.svd(C); rank=np.sum(sv>1e-9*sv[0]) if len(sv) else 0
    Kb=vt[rank:].T                                   # columns span ker(C)
    Gp=Kb.T@Gw@Kb; Mp=Kb.T@Mm@Kb
    ev,ec=np.linalg.eigh(np.linalg.solve(Mp,Gp) if False else np.linalg.pinv(Mp)@Gp)
    order=np.argsort(ev.real); ev=ev.real[order]; ec=ec.real[:,order]
    cmin=Kb@ec[:,0]; cmin/=np.linalg.norm(cmin)
    return ev[0], cmin

print("Constrained Weil minimizer over prolate modes {h_0,h_2,...}, s.t. f(0)=fhat(0)=0. Does it -> k_lambda (h_0,h_4)?")
for lmbda in [1.6, 2.0]:
    Gw,Mm,h0,intg,cols=analyze(lmbda)
    lab=[f"h{2*j}" for j in range(len(cols))]
    for cond,C in [("fhat(0)=0 only",np.array([intg])),
                   ("f(0)=fhat(0)=0",np.array([h0,intg]))]:
        val,cmin=constrained_min(Gw,Mm,C)
        # normalize sign, show dominant components
        cmin=cmin/np.max(np.abs(cmin))
        comp=", ".join(f"{lab[i]}:{cmin[i]:+.2f}" for i in range(len(cols)))
        print(f"  lambda={lmbda} [{cond:>16}] minQW={val:>11.3e}  dir=({comp})",flush=True)
print("\nREAD: under f(0)=fhat(0)=0, is the minimizer DOMINATED by h0,h4 (k_lambda)? If yes = the 2 conditions SELECT")
print("k_lambda as the near-null Weil direction = the analytic core of Connes Step 2 (well-conditioned, floor-immune).")
