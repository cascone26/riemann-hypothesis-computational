#!/usr/bin/env python3
"""
g_a_compact.py — LEDGER §110: Suzuki §8 COMPACT reformulation of λ_a, to break the large-a wall (§42/§109).
Instead of the sine-Galerkin of the unbounded, pole-SINGULAR A_a, use Suzuki §8:
  λ_a = inf_{v∈H^1_0(-a,a)} Q_W(v)/||v||²,  Q_W(v)=∫∫ g(x-y) v'(x)v'(y).
Substitute u = v' (mean-zero, since v(±a)=0): Q_W = <G u,u>, G kernel g(x-y); and ||v||² = <N u,u>,
where N = (DD*)^{-1} = inverse Neumann Laplacian, EXPLICIT continuous kernel (§8.2):
  N(x,y) = (x²+y²)/(4a) − |x−y|/2 + a/6.
⟹ λ_a = smallest generalized eigenvalue of  G u = λ N u  over mean-zero u.
Both G and N are integral operators with CONTINUOUS kernels (Nyström discretization) — much better conditioned than
the derivative-based Galerkin of the singular A_a. Test: does this reach large a cleanly where A_a's λ_a went
spuriously negative (§42) / floored? Also cross-check vs the mpmath A_a values (§109): a=0.6→+2.3e-5 etc.
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def lambda_a_compact(a, M=1200):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    # G kernel: g(x-y) (Toeplitz), g = full screw
    dif=np.arange(-M,M+1)*h
    gv=np.array([screw_g(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    # N kernel (inverse Neumann Laplacian), explicit continuous
    X=x[:,None]; Y=x[None,:]
    N=(X**2+Y**2)/(4*a) - np.abs(X-Y)/2 + a/6
    # quadrature-weighted operator matrices (symmetric): A_ij = sqrt(w_i) K(x_i,x_j) sqrt(w_j)
    sw=np.sqrt(w)
    Gm=(sw[:,None]*G*sw[None,:]); Gm=(Gm+Gm.T)/2
    Nm=(sw[:,None]*N*sw[None,:]); Nm=(Nm+Nm.T)/2
    # enforce mean-zero u: project onto complement of the (weighted) constant vector
    c=sw.copy(); c/=np.linalg.norm(c)                 # weighted constant direction
    P=np.eye(M+1)-np.outer(c,c)
    Gp=P@Gm@P; Np=P@Nm@P
    Gp=(Gp+Gp.T)/2; Np=(Np+Np.T)/2
    # generalized eig on the mean-zero subspace: use pseudo-inverse of Np (it's singular in the c-direction)
    # solve Np^{-1/2} Gp Np^{-1/2} eigenproblem on the range of P
    ev_N,V_N=np.linalg.eigh(Np)
    keep=ev_N>1e-12*ev_N.max()
    S=V_N[:,keep]/np.sqrt(ev_N[keep])                 # Np^{-1/2} on its range
    Ared=S.T@Gp@S; Ared=(Ared+Ared.T)/2
    lam=np.linalg.eigvalsh(Ared)
    return lam[0], lam[1]-lam[0]

print("Suzuki §8 COMPACT reformulation: λ_a via  G u = λ N u  (continuous kernels). Test large-a conditioning.")
print(f"{'a':>5} {'λ_a(compact)':>16} {'gap':>12}   (cf mpmath A_a §109: 0.6→+2.3e-5,1.0→+1.6e-5,2.0→+1.3e-5)")
for a in [0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.5, 3.0]:
    try:
        la,gap=lambda_a_compact(a)
        flag='' if la>-1e-9 else '  <-- NEGATIVE'
        print(f"{a:>5.1f} {la:>16.6e} {gap:>12.4e}{flag}",flush=True)
    except Exception as e:
        print(f"{a:>5.1f}  ERROR: {e}",flush=True)
print("\nREAD: does λ_a stay POSITIVE + clean into large a (2.5,3.0) where the sine-Galerkin A_a failed (§42)?")
print("If the compact G_a reformulation is better-conditioned, it reaches the zeta regime the singular A_a couldn't.")
