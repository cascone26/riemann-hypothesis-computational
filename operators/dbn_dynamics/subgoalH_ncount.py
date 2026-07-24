#!/usr/bin/env python3
"""
subgoalH_ncount.py — LEDGER §134 (Sub-goal H, the reweighting/reformulation test): measure the GROWTH of the near-null
core n_neg(a) of the Weil form, and whether a REFORMULATION (compact G_a §110, prolate-conditioned) bounds/slows it.
CCM claim: n_neg ~ 1.8 log c grows with cutoff c; Sub-goal H = a reweighting keeps it uniformly finite + margin held.
Sylvester: the COUNT is congruence-invariant, so the real question is the GROWTH LAW of the near-null cluster and its
conditioning (the margin gap) under different formulations. We compute the full low spectrum vs a and count the cluster.
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def weil_spectrum_compact(a, M=1400):
    """low spectrum of the Weil form via the compact G_a formulation (§110)."""
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h; gv=np.array([screw_g(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    X=x[:,None]; Y=x[None,:]; N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w); Gm=sw[:,None]*G*sw[None,:]; Gm=(Gm+Gm.T)/2; Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c); Gp=P@Gm@P; Np=P@Nm@P
    ev_N,V_N=np.linalg.eigh((Np+Np.T)/2); keep=ev_N>1e-12*ev_N.max()
    S=V_N[:,keep]/np.sqrt(ev_N[keep]); Ar=S.T@((Gp+Gp.T)/2)@S
    return np.linalg.eigvalsh((Ar+Ar.T)/2)

print("Sub-goal H: near-null CORE growth n_neg(a) of the Weil form (compact G_a formulation). Cutoff c ~ e^{2a}.")
print(f"{'a':>5} {'log c=2a':>9} {'λ_min':>11} {'#(<1e-3)':>9} {'#(<1e-2)':>9} {'#(<0.1)':>8} {'bulk starts~':>12}")
rows=[]
for a in [0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6,4.0]:
    ev=weil_spectrum_compact(a)
    lo=ev[ev<0.5]                                       # near-null + intermediate, below the O(1) bulk
    n3=int(np.sum(ev<1e-3)); n2=int(np.sum(ev<1e-2)); n1=int(np.sum(ev<0.1))
    # bulk start = first eigenvalue that's O(1) (jump)
    bulk = ev[np.argmax(ev>0.3)] if np.any(ev>0.3) else np.nan
    rows.append((2*a, n3, n2, n1))
    print(f"{a:>5.1f} {2*a:>9.2f} {ev[0]:>11.3e} {n3:>9} {n2:>9} {n1:>8} {bulk:>12.3f}",flush=True)
# growth law: fit n(<1e-2) vs log c = 2a
lc=np.array([r[0] for r in rows]); n2v=np.array([r[2] for r in rows])
A=np.polyfit(lc, n2v, 1)
print(f"\n GROWTH: #(<1e-2) ≈ {A[0]:.2f}·(log c) + {A[1]:.2f}   [CCM predict slope ~1.8]")
print(" READ: does the near-null core GROW ~log c (⟹ finite-core reduction needs the growing piece handled) — and does")
print(" the compact reformulation's slope beat the naive ~1.8? A bounded/slower core = the Sub-goal H signal.")
