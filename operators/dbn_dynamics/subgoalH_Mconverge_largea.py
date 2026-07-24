#!/usr/bin/env python3
"""
subgoalH_Mconverge_largea.py — LEDGER §139 (Sub-goal H, resolving the §138 large-a artifact): §138's λ_min sweep
looked like it DECAYED past a=3 (9e-7 → 4.8e-7 at a=5), BUT its own conv% column flagged the bottom mode as
UN-resolved there (M=1800 vs M=2600 disagreed 2-4x at a≥3.5). §137 only licensed grid-trust up to a=3. So the
"decay" is almost certainly under-resolution, not marginality. RESOLVE IT: fix the worst-case large a and push M up
until λ1 stops moving. If λ1 converges back to ~9e-7 (the small-a floor) ⟹ λ_min is FLAT/uniform (RH-favorable). If
it converges to something clearly smaller ⟹ real decay (marginal). tr(N)~a² (§138) already ruled out a hidden e^{a}.
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def lam_low(a, M, k=4):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h; gv=np.array([screw_g(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    X=x[:,None]; Y=x[None,:]; N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w); Gm=sw[:,None]*G*sw[None,:]; Gm=(Gm+Gm.T)/2; Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c); Gp=P@Gm@P; Np=P@Nm@P
    ev_N,V_N=np.linalg.eigh((Np+Np.T)/2); keep=ev_N>1e-12*ev_N.max()
    S=V_N[:,keep]/np.sqrt(ev_N[keep]); Ar=S.T@((Gp+Gp.T)/2)@S
    return np.linalg.eigvalsh((Ar+Ar.T)/2)[:k]

print("Sub-goal H §139: M-convergence of λ1 at fixed LARGE a (resolve §138's suspected under-resolution artifact).")
print("If λ1 climbs back to the small-a floor ~9e-7 as M↑ ⟹ λ_min FLAT/uniform (RH-favorable). Settles lower ⟹ real decay.\n")
Ms=[1800,2600,3600,5000,6500]
for a in [4.0, 4.5, 5.0]:
    print(f"=== a={a} (log c={2*a:.1f}) ===")
    print(f"  {'M':>6} {'λ1':>12} {'λ2':>12} {'λ3':>12} {'Δλ1 vs prev':>13}")
    prev=None
    for M in Ms:
        lv=lam_low(a,M)
        d = '' if prev is None else f"{(lv[0]-prev)/prev*100:+.1f}%"
        print(f"  {M:>6} {lv[0]:>12.4e} {lv[1]:>12.4e} {lv[2]:>12.4e} {d:>13}",flush=True)
        prev=lv[0]
    print()
print("READ: watch the Δλ1 column shrink toward 0 (convergence). The converged λ1 is the TRUE λ_min at that a.")
print("Compare to the small-a floor 9.2e-7: match ⟹ uniform positivity confirmed out to a=5 ⟹ Sub-goal H favorable HERE.")
