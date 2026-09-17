"""
Test the explicit near-optimal function psi(x)=cosh(x/2)-cosh(a/2) (admissible: vanishes at +-a, cosh-like,
= the negative-square direction phi_1). Compute R(a,psi)=int int g(x-y)psi'(x)psi'(y)/||psi||^2, psi'=1/2 sinh(x/2).
Compare to lambda_a. If R(a,psi) ~ lambda_a and >=0, cosh is the explicit near-optimal test function.
"""
import numpy as np
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto, lambda_a

def R_psi(a, M):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gv=np.array([screw_g(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1); G=gv[(idx[:,None]-idx[None,:])+M]
    psi=np.cosh(x/2)-np.cosh(a/2)
    dpsi=0.5*np.sinh(x/2)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dpw=dpsi*w
    num=dpw@G@dpw
    nrm=np.sum(psi*psi*w)         # ||psi||^2 L2
    return num/nrm

if __name__=="__main__":
    print("Explicit test fn psi=cosh(x/2)-cosh(a/2): R(a,psi) vs lambda_a")
    print(f"{'a':>5} {'R(a,psi)':>13} {'lambda_a':>13} {'R>=0?':>6}")
    for a in [0.5,1.0,2.0,3.0]:
        M=min(3000,max(1500,int(400/a)))
        Rp=R_psi(a,M); la=lambda_a(a,M=M,K=80)
        print(f"{a:5.2f} {Rp:13.6e} {la:13.6e} {'OK' if Rp>=0 else 'NEG':>6}",flush=True)
    print("\nR(a,psi) is an UPPER bound on lambda_a (single test fn). If ~lambda_a and >=0, cosh is near-optimal.")
