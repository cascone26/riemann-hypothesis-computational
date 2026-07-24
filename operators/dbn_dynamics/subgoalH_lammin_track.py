#!/usr/bin/env python3
"""
subgoalH_lammin_track.py — LEDGER §138 (Sub-goal H, the decisive measurement): now that §137 PROVED the near-null
modes are grid-independent (best-overlap→1 for the lowest ~50 modes across 3x refinement at all a≤3), the eigenvalues
are trustworthy — so track λ_min(a) and the bottom ladder OUT to large a and ask the uniform-positivity question.

§137 surprise: the low eigenvalues were nearly a-INDEPENDENT over a∈[1.5,3.0] (λ1≈9.2e-7 at all four). If that holds
out to a=5, the compact form is uniformly positive with a scale-universal near-null ladder ⟹ Sub-goal H favorable.
But we must RULE OUT a hidden a-scale in the Neumann-mass normalization N before calling it "universal":
  - report λ1,λ2,λ3 with per-point grid convergence (M=1800 vs 2600) so we trust each value,
  - report the RAW G-Rayleigh of the bottom mode (physical scale), plus tr(N) vs a (hidden e^{a} factor?).
DECISIVE: λ_min(a) flat as a→∞ ⟹ uniform positivity (RH-favorable, tractable core). λ_min(a)→0 ⟹ marginal (the
reweighting must actively hold the margin — the hard CCM case).
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def weil_full(a, M):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h; gv=np.array([screw_g(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    X=x[:,None]; Y=x[None,:]; N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w); Gm=sw[:,None]*G*sw[None,:]; Gm=(Gm+Gm.T)/2; Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c); Gp=P@Gm@P; Np=P@Nm@P
    ev_N,V_N=np.linalg.eigh((Np+Np.T)/2); keep=ev_N>1e-12*ev_N.max()
    S=V_N[:,keep]/np.sqrt(ev_N[keep]); Ar=S.T@((Gp+Gp.T)/2)@S
    wv,U=np.linalg.eigh((Ar+Ar.T)/2)
    v0=S@U[:,0]; gray=float(v0@Gp@v0); nray=float(v0@Np@v0)   # bottom-mode raw G- and N-Rayleigh (nray≈1)
    N_scale=float(np.trace(Np))
    return wv, gray, nray, N_scale

print("Sub-goal H §138: λ_min(a) and near-null ladder OUT to large a — is the compact Weil form UNIFORMLY positive?")
print("(§137 proved these eigenvalues are grid-trustworthy. Question: λ_min flat ⟹ RH-favorable; →0 ⟹ marginal.)\n")
print(f"{'a':>4} {'logc':>5} {'λ1(M1800)':>11} {'λ1(M2600)':>11} {'conv%':>6} {'λ2':>10} {'λ3':>10} {'λ2/λ1':>7} {'tr(N)':>10}")
lam1=[]
for a in [1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]:
    w1,_,_,_ = weil_full(a,1800)
    w2,gray,nray,Nsc = weil_full(a,2600)
    conv = 100*(1-abs(w1[0]-w2[0])/abs(w2[0])) if w2[0]!=0 else 0
    lam1.append((2*a,w2[0]))
    print(f"{a:>4.1f} {2*a:>5.1f} {w1[0]:>11.3e} {w2[0]:>11.3e} {conv:>6.1f} {w2[1]:>10.3e} {w2[2]:>10.3e} {w2[1]/w2[0]:>7.2f} {Nsc:>10.3e}",flush=True)
print("\nλ_min(a) trend:")
lc=np.array([r[0] for r in lam1]); l1=np.array([r[1] for r in lam1])
for (c,v) in lam1:
    print(f"   logc={c:.1f}: λ1={v:.4e}")
sl_log=np.polyfit(lc, np.log(l1), 1)[0]
sl_lin=np.polyfit(lc, l1, 1)[0]
print(f"\n   d(log λ1)/d(log c) = {sl_log:.3f}   (≈0 ⟹ λ_min FLAT = uniform positivity; negative ⟹ decaying = marginal)")
print(f"   d(λ1)/d(log c)     = {sl_lin:.2e}")
print("READ: if λ1 stays ~1e-6 flat out to a=5 AND tr(N) isn't secretly absorbing an e^{a} ⟹ the compact Weil form is")
print("UNIFORMLY POSITIVE with a scale-universal near-null floor — Sub-goal H's 'uniform finiteness+margin' holds HERE,")
print("moving the open question to whether this compact positivity transfers back to the true (non-compact) Weil form.")
