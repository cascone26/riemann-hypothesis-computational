#!/usr/bin/env python3
"""
subgoalH_ghost_discriminator.py — LEDGER §140 (Sub-goal H, the crux discriminator): §139 showed λ_min at a=4 DROPPING
with grid M (1.03e-6→1.85e-7, not converging). Two mutually-exclusive explanations, opposite meaning:
  (A) REAL: at a>3 the true near-null core (§134 doubly-exp) reasserts — finer grids resolve genuinely smaller Weil modes.
  (B) GHOST: the N-whitening (S = V_N/sqrt(ev_N), keep>1e-12·max) manufactures spurious tiny modes as N gets
      ill-conditioned at large M — pure numerical artifact of the compact reformulation, not Weil structure.
DISCRIMINATE with two independent cross-checks at the clearest problem point (a=4):
  1. KEEP-THRESHOLD SWEEP: recompute λ1 with keep ∈ {1e-13,1e-11,1e-9,1e-7}. Ghost ⟹ λ1 shrinks as keep loosens
     (more near-null N-directions admitted). Real ⟹ λ1 stable across thresholds.
  2. DIRECT GENERALIZED SOLVER: scipy.linalg.eigh(Gp, Np) (symmetric-definite generalized, no manual whitening) vs the
     manual-whitening λ1. Agree ⟹ real. Disagree (direct gives larger/stable λ1) ⟹ whitening artifact.
  3. cond(N_kept) reported — if it explodes (>1e12) at large M, whitening is untrustworthy there by definition.
"""
import numpy as np, sys
from scipy.linalg import eigh as geigh
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def build_GN(a, M):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h; gv=np.array([screw_g(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    X=x[:,None]; Y=x[None,:]; N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w); Gm=sw[:,None]*G*sw[None,:]; Gm=(Gm+Gm.T)/2; Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c)
    Gp=P@Gm@P; Np=P@Nm@P; Gp=(Gp+Gp.T)/2; Np=(Np+Np.T)/2
    return Gp, Np

def lam_whiten(Gp, Np, keep_rel):
    ev_N,V_N=np.linalg.eigh(Np); mx=ev_N.max(); keep=ev_N>keep_rel*mx
    cond = mx/ev_N[keep].min()
    S=V_N[:,keep]/np.sqrt(ev_N[keep]); Ar=S.T@Gp@S
    return np.linalg.eigvalsh((Ar+Ar.T)/2)[0], int(keep.sum()), cond

def lam_direct(Gp, Np):
    # symmetric-definite generalized eigenproblem on the range of Np (drop the 1 null direction from the P-projection)
    ev_N,V_N=np.linalg.eigh(Np); keep=ev_N>1e-10*ev_N.max()
    Q=V_N[:,keep]                                   # restrict to Np's positive range (kills the projected-out null dir)
    Gr=Q.T@Gp@Q; Nr=Q.T@Np@Q; Gr=(Gr+Gr.T)/2; Nr=(Nr+Nr.T)/2
    wv=geigh(Gr, Nr, eigvals_only=True)             # direct generalized solver, no sqrt-whitening
    return wv[0]

print("Sub-goal H §140: is the large-a λ_min drop REAL (doubly-exp core) or a GHOST (N-whitening conditioning)?")
print("Discriminators at a=4.0: keep-threshold sweep + direct generalized solver + cond(N).\n")
for a,M in [(3.0,2600),(4.0,2600),(4.0,3600),(4.0,5000),(4.5,3600)]:
    Gp,Np = build_GN(a,M)
    print(f"--- a={a}, M={M} ---")
    print(f"  {'keep_rel':>10} {'λ1(whiten)':>13} {'#kept':>7} {'cond(N)':>11}")
    for kr in [1e-13,1e-11,1e-9,1e-7]:
        l1,nk,cd = lam_whiten(Gp,Np,kr)
        print(f"  {kr:>10.0e} {l1:>13.4e} {nk:>7} {cd:>11.2e}",flush=True)
    ld = lam_direct(Gp,Np)
    print(f"  {'DIRECT eigh(G,N)':>10}: λ1 = {ld:>.4e}",flush=True)
    print()
print("READ: if λ1(whiten) is STABLE across keep_rel AND ≈ λ1(direct) ⟹ the drop is REAL (true small Weil modes at a>3).")
print("If λ1 SHRINKS as keep_rel loosens and/or DISAGREES with direct ⟹ GHOST (whitening artifact); the real λ_min is the")
print("stable/direct value. Either way this tells us whether §138's flat-9e-7-floor breaks for real past a=3, or just numerically.")
