# Connes program crux (LEDGER §94/96), STAGE 3 (the OPEN estimate): compare theta_x = Suzuki-screw minimal
# eigenvector to k_lambda = E(prolate combo). Connes Step 2: theta_x ~ k_lambda. We already know even-simple
# holds (S48/49); the PRIZE here is the STRUCTURE of the correction theta_x - k_lambda -- what an analytic proof needs.
# a = log(lambda). Trust a<=~0.9 (suzuki flags a>1 as discretization artifact). Compare unit-normalized, sign-aligned,
# on (-a,a). Test BOTH v and v'=f against k_lambda (convention: Suzuki form uses v', norm uses v).
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto
from numpy.polynomial.legendre import legval

def theta_suzuki(a, M=3000, K=60):
    """minimal eigenVECTOR of the Suzuki screw Weil form on (-a,a). Returns v(x) and v'(x) on grid x."""
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gvals=np.array([screw_g(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1); G=gvals[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1)
    arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg)                    # basis derivatives
    Vv=np.sin(arg)                                     # basis functions
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]
    Qmat=Wp.T@G@Wp; Qmat=(Qmat+Qmat.T)/2
    ev,ec=np.linalg.eigh(Qmat/a)
    c=ec[:,0]                                          # minimal eigenvector coefficients
    v=Vv@c; vp=Vp@c                                    # v(x), v'(x)
    return x, v, vp, ev[0]

# ---- prolate k_lambda (from k_lambda_construct) ----
def y2_matrix(Nleg):
    Y=np.zeros((Nleg+2,Nleg+2))
    for n in range(Nleg+1):
        if n+1<=Nleg+1: Y[n+1,n]=(n+1)/(2*n+1)*np.sqrt((2*n+1)/(2*n+3))
        if n-1>=0:      Y[n-1,n]=n/(2*n+1)*np.sqrt((2*n+1)/(2*n-1))
    return (Y[:Nleg,:Nleg])@(Y[:Nleg,:Nleg])
def k_lambda_on(s_grid, lmbda):
    Nleg=max(200,int(6*2*np.pi*lmbda**2)); c=2*np.pi*lmbda**2; n=np.arange(Nleg)
    A=np.diag((n*(n+1)).astype(float))+c**2*y2_matrix(Nleg); A=(A+A.T)/2
    w,V=np.linalg.eigh(A)
    par=np.array([(-1)**int(np.argmax(np.abs(V[:,k]))) for k in range(V.shape[1])])
    ei=np.where(par>0)[0]; i0,i4=ei[0],ei[2]
    a0,b0=1.0,(-V[0,i0]/V[0,i4] if abs(V[0,i4])>1e-30 else 0.0)
    def coeff(col): return np.array([col[m]*np.sqrt((2*m+1)/2) for m in range(Nleg)])
    c0,c4=coeff(V[:,i0]),coeff(V[:,i4])
    def h_lam(xv):
        y=xv/lmbda; out=np.zeros_like(xv,float); m=np.abs(y)<=1
        out[m]=a0*legval(y[m],c0)+b0*legval(y[m],c4); return out
    u=np.exp(s_grid); nmax=int(np.floor(lmbda**2))+1
    ssum=np.zeros_like(u)
    for nn in range(1,nmax): ssum+=h_lam(nn*u)
    return np.sqrt(u)*ssum

def unit(z):
    z=z/np.sqrt(np.trapezoid(z*z)); return z
def align_dist(A,B):
    A,B=unit(A.copy()),unit(B.copy())
    if np.trapezoid(A*B)<0: B=-B
    return np.sqrt(np.trapezoid((A-B)**2)), np.trapezoid(A*B)   # L2 dist, overlap

print("STAGE 3: Suzuki minimal eigenvector theta vs Connes k_lambda, a=log(lambda). (unit-norm, sign-aligned on (-a,a))\n")
print(f"{'lambda':>7} {'a':>6} {'lam_a':>10} {'parity(v)':>10} {'dist(v,k)':>10} {'ov(v,k)':>8} {'dist(vp,k)':>11} {'ov(vp,k)':>9}")
for lmbda in [1.8, 2.0, 2.2, 2.45]:
    a=np.log(lmbda)
    x,v,vp,la=theta_suzuki(a)
    k=k_lambda_on(x,lmbda)
    # parity of v about x=0
    parity='even' if np.trapezoid(v*v[::-1])>0.9*np.trapezoid(v*v) else ('odd' if np.trapezoid(v*v[::-1])<-0.9*np.trapezoid(v*v) else 'mixed')
    dvk,ovk=align_dist(v,k)
    dpk,opk=align_dist(vp,k)
    print(f"{lmbda:>7.2f} {a:>6.3f} {la:>10.3e} {parity:>10} {dvk:>10.4f} {ovk:>8.4f} {dpk:>11.4f} {opk:>9.4f}",flush=True)
print("\nREAD: which of v / v' matches k_lambda (overlap->1, dist->0)? That fixes the convention. Then the RESIDUAL")
print("structure (dist not ->0) = the correction theta-k_lambda whose FORM the analytic Step-2 proof must control.")
