#!/usr/bin/env python3
"""
subgoalH_jacobi_structural.py — LEDGER §147 (Lens 1 ∩ Lens 7 ∩ §119 convergence): the compact Weil form's bottom is a
JACOBI operator (§146). Extract its three-term-recurrence coefficients and ask THE structural question: is its
positivity λ_min>0 certified by STRUCTURE (diagonal dominance / oscillation / total positivity) rather than by computing
all eigenvalues? A structural certificate that persists as a→∞ would be a positivity PROOF mechanism (→ RH); if instead
the structural margin collapses exactly like λ_min (doubly-exp), the Jacobi form faithfully carries the knife-edge.

Build the compact Weil operator in the sine VALUE basis: A_kl=∫∫g(x-y)v_k v_l, B_kl=∫∫N(x,y)v_k v_l (Neumann mass); the
symmetric operator C=B^{-1/2}A B^{-1/2} has the λ_k spectrum. Test:
  (1) band structure of C in the sine basis (confirm §144 tridiagonal-on-parity = Jacobi).
  (2) restrict to a parity class → tridiagonal J; extract diag α_k, offdiag β_k.
  (3) GERSHGORIN / diagonal-dominance margin: min_k(α_k − |β_{k-1}| − |β_k|). >0 ⟹ structural positivity certificate.
  (4) compare that structural margin to the true λ_min, across a — does structure EXPLAIN the positivity, and how does
      the margin scale?
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def sine_value_C(a, gfun, K=24, M=2400):
    """CORRECT Weil operator in the sine basis: uses DERIVATIVES v_k' (the mean-zero structure that makes the Weil form
    positive; §144). Q_kl=∫∫g(x-y)v_k'(x)v_l'(y), operator C=Q/a (mass ‖v‖²=aΣc²). Eigenvalues = λ_k>0."""
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h; gv=np.array([gfun(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); Vp=(kk*np.pi/(2*a))*np.cos(np.pi*np.outer((x+a),kk)/(2*a))   # v_k'(x) (M+1,K)
    Wp=Vp*w[:,None]
    Q=Wp.T@G@Wp; Q=(Q+Q.T)/2
    return Q/a

def band_fracs(C):
    K=C.shape[0]; k=np.arange(K); fro=np.sqrt(np.sum(C**2))
    return {b: round(np.sqrt(np.sum((C*(np.abs(np.subtract.outer(k,k))<=b))**2))/fro,3) for b in [0,1,2,3]}

print("Sub-goal H §147: is the bottom JACOBI operator's positivity STRUCTURAL? (diag-dominance vs true λ_min)\n")
for a in [1.5,2.0,2.5,3.0]:
    C=sine_value_C(a, screw_g)
    ev=np.linalg.eigvalsh(C)
    bf=band_fracs(C)
    # even sublattice (indices 0,2,4,... in 0-based = sine modes k=1,3,5..) — pick the parity block that holds the ground mode
    for par,name in [(0,"even-idx"),(1,"odd-idx")]:
        sub=C[par::2][:,par::2]; K2=sub.shape[0]
        al=np.diag(sub); be=np.diag(sub,1)
        # diagonal-dominance (Gershgorin) lower bound on λ_min of this tridiag-ish block
        radius=np.abs(np.concatenate([[0],be]))+np.abs(np.concatenate([be,[0]]))
        gersh=np.min(al-radius)
        subev=np.linalg.eigvalsh(sub)
        # tridiagonality of the block: fraction of energy on |i-j|<=1
        kk=np.arange(K2); frob=np.sqrt(np.sum(sub**2)); tri=np.sqrt(np.sum((sub*(np.abs(np.subtract.outer(kk,kk))<=1))**2))/frob
        print(f" a={a} [{name}] λ_min(block)={subev[0]:+.3e} Gershgorin={gersh:+.3e} tridiag-frac={tri:.3f} "
              f"α_1..3={np.round(al[:3],3)} β_1..3={np.round(be[:3],3)}",flush=True)
    print(f"     full C: λ_min={ev[0]:+.3e}  band-frac|k-l|≤{{0,1,2,3}}={bf}",flush=True)
    print()
print("READ: tridiag-frac→1 ⟹ genuine Jacobi. Gershgorin>0 ⟹ STRUCTURAL positivity certificate (diag dominance proves")
print("λ_min>0 without eigen-solve). If Gershgorin>0 and tracks λ_min as a grows ⟹ a candidate structural positivity")
print("mechanism (Lens 7). If Gershgorin<0 while λ_min>0 ⟹ positivity is NON-structural (cancellation), no cheap cert.")
