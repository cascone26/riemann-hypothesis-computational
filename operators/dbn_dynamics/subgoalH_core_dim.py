#!/usr/bin/env python3
"""
subgoalH_core_dim.py — LEDGER §136 (Sub-goal H, resolving the §135 ambiguity): does the near-null CORE of the Weil
form have a STABLE FINITE DIMENSION, or does it genuinely proliferate?

§135 counted eigenvalues below a threshold and saw the count grow with grid resolution M, then RETRACTED the "bounded
core" as a grid artifact. But eigenvalue-COUNTING cannot tell "finite core + grid noise near zero" from "genuinely
infinite core" — both make the count grow with M. The decisive test is on the EIGENVECTORS:

  Compute the near-null mode SHAPES at several M, interpolate to a common reference grid, and measure the principal
  angles between the lowest-k subspace at M and at M+ΔM. If the shapes STABILIZE (angles→0) for the lowest k modes and
  only THEN blow up, the true core is exactly k-dimensional and FINITE (Sub-goal H tractable). If every mode churns
  with M, the core is genuinely infinite.

Relative subspace overlaps are numerically stable even when the absolute eigenvalues are ~1e-295 (doubly-exp), so this
sidesteps the precision wall that made §135's count approach hopeless. Local float64, fast.
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def weil_modes_compact(a, M, ncore=40):
    """low spectrum AND mode shapes of the compact Weil form G_a. Returns (evals[:ncore], funcs[:,:ncore] on x-grid, x)."""
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h; gv=np.array([screw_g(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    X=x[:,None]; Y=x[None,:]; N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w); Gm=sw[:,None]*G*sw[None,:]; Gm=(Gm+Gm.T)/2; Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c); Gp=P@Gm@P; Np=P@Nm@P
    ev_N,V_N=np.linalg.eigh((Np+Np.T)/2); keep=ev_N>1e-12*ev_N.max()
    S=V_N[:,keep]/np.sqrt(ev_N[keep]); Ar=S.T@((Gp+Gp.T)/2)@S
    wv,U=np.linalg.eigh((Ar+Ar.T)/2)                       # ascending eigenvalues + eigvecs of the whitened form
    funcs = S@U                                            # (M+1, k) generalized eigenvectors as mode shapes on x-grid
    # de-whiten to actual x-space amplitude: v_whitened lives in sqrt(w)-space; divide out sw to get the raw function
    funcs = funcs/np.where(sw[:,None]>0,sw[:,None],1.0)
    return wv[:ncore], funcs[:,:ncore], x

def subspace_angles(F1,x1,F2,x2,k,xref):
    """principal angles between span(lowest-k modes at M1) and span(lowest-k modes at M2), on a common grid xref."""
    # interpolate each mode to the common reference grid, L2-orthonormalize each subspace, then SVD the cross-Gram
    def interp_ortho(F,x):
        A=np.array([np.interp(xref,x,F[:,j]) for j in range(k)]).T   # (len(xref), k)
        Q,_=np.linalg.qr(A)                                          # orthonormal columns (L2 on uniform xref ~ up to h)
        return Q
    Q1=interp_ortho(F1,x1); Q2=interp_ortho(F2,x2)
    s=np.linalg.svd(Q1.T@Q2, compute_uv=False)
    s=np.clip(s,-1,1)
    return np.degrees(np.arccos(s))     # principal angles in degrees; small = subspaces agree

print("Sub-goal H §136: is the near-null CORE a STABLE FINITE-dim space, or does it proliferate? (eigenVECTOR test)")
print("Method: principal angles between lowest-k mode subspace at grid M vs M+ΔM. Angles→0 for k≤k* then blow up ⟹")
print("        the true core is exactly k*-dimensional and FINITE. Every mode churning ⟹ genuinely infinite core.\n")

for a in [1.5, 2.0, 2.5]:
    Ms=[600,900,1200,1600]
    xref=np.linspace(-a,a,300)
    data=[weil_modes_compact(a,M) for M in Ms]
    print(f"=== a={a}  (log c = {2*a:.1f}) ===")
    print(f"  lowest 8 eigenvalues at M={Ms[-1]}: "+", ".join(f"{v:.2e}" for v in data[-1][0][:8]))
    # for growing k, report the MAX principal angle between the FINEST two grids (most converged comparison)
    (w1,F1,x1),(w2,F2,x2)=data[-2],data[-1]
    print(f"  {'k':>3} {'max principal angle (deg), M=1200 vs 1600':>44}   {'verdict':>10}")
    for k in [1,2,3,4,5,6,8,10,14,20]:
        ang=subspace_angles(F1,x1,F2,x2,k,xref)
        verdict='STABLE' if ang.max()<15 else ('drifting' if ang.max()<45 else 'CHURN')
        print(f"  {k:>3} {ang.max():>44.2f}   {verdict:>10}",flush=True)
    print()
print("READ: the largest k whose subspace stays STABLE (max angle small) across refinement = the true finite core dim.")
print("If that k* is fixed and small while eigenvalue COUNTS grow ⟹ §135's growth WAS grid noise on a finite core ⟹")
print("Sub-goal H reweighting only needs to handle a fixed finite-dim space (tractable). If k*=0/churn ⟹ genuinely infinite.")
