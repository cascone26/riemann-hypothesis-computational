#!/usr/bin/env python3
"""
subgoalH_core_dim_v2.py — LEDGER §137 (Sub-goal H, decisive version): count the STABLE near-null core dimension in a
way that is IMMUNE to eigenvalue-ordering swaps — the confound that made §136's fixed-integer-k comparison
non-monotone (a=2.0: k=10 churned, k=14 re-stabilized ⟹ near-degenerate modes swapping order between grids).

FIX: don't compare "lowest-k subspace at M1 vs M2" by index. Instead, for EACH near-null mode at grid M1, find its
BEST overlap (max |cos| of principal angle) against the ENTIRE near-null space at the finer grid M2. A mode that is a
REAL, grid-independent core mode has overlap → 1; a grid-noise mode has overlap < 1 and no stable partner. The count
  d*(M) = #{ modes at M with best-overlap-to-finer-grid > 0.97 }
is reordering-proof. If d*(M) CONVERGES to a fixed value as M grows, that value is the true finite core dimension at
this a. Then we watch d* vs a: BOUNDED ⟹ Sub-goal H reweighting handles a fixed finite space (tractable);
GROWING without bound ⟹ genuinely infinite core.
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def weil_modes(a, M, ncore=60):
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
    funcs=(S@U); funcs=funcs/np.where(sw[:,None]>0,sw[:,None],1.0)
    return wv[:ncore], funcs[:,:ncore], x

def best_overlaps(F1,x1,F2,x2,nq, nspace, xref):
    """for each of the lowest nq modes at grid-1, its best |cos| overlap against the lowest-nspace span at grid-2."""
    A1=np.array([np.interp(xref,x1,F1[:,j]) for j in range(nq)]).T        # query modes (len xref, nq)
    A2=np.array([np.interp(xref,x2,F2[:,j]) for j in range(nspace)]).T    # reference span (len xref, nspace)
    Q2,_=np.linalg.qr(A2)                                                 # orthonormal basis of finer-grid low space
    # normalize each query mode, project onto Q2, overlap = ||projection|| (=cos of angle to the subspace)
    ov=[]
    for j in range(nq):
        v=A1[:,j]; v=v/np.linalg.norm(v)
        ov.append(np.linalg.norm(Q2.T@v))
    return np.array(ov)

print("Sub-goal H §137 (decisive): reordering-proof stable-core dimension d* = #{modes with best-overlap-to-finer-grid>0.97}.")
print("Converges in M ⟹ that's the true finite core dim at this a. Then: is d* BOUNDED or GROWING in a?\n")

results={}
for a in [1.5, 2.0, 2.5, 3.0]:
    xref=np.linspace(-a,a,400)
    Ms=[800,1200,1800,2600]
    modes=[weil_modes(a,M) for M in Ms]
    print(f"=== a={a}  (log c={2*a:.1f}) ===")
    print(f"  lowest 6 λ (M={Ms[-1]}): "+", ".join(f"{v:.2e}" for v in modes[-1][0][:6]))
    dstars=[]
    for i in range(len(Ms)-1):
        (w1,F1,x1)=modes[i]; (w2,F2,x2)=modes[i+1]
        ov=best_overlaps(F1,x1,F2,x2, nq=50, nspace=55, xref=xref)
        for thr in [0.99,0.97]:
            pass
        d97=int(np.sum(ov>0.97)); d99=int(np.sum(ov>0.99))
        dstars.append(d97)
        # find first mode index that FALLS below 0.97 (the edge of the stable core)
        below=np.where(ov<=0.97)[0]; edge=below[0] if len(below) else 50
        print(f"  M {Ms[i]:>4}->{Ms[i+1]:>4}:  d*(>.99)={d99:>3}  d*(>.97)={d97:>3}  first-unstable-mode@={edge:>3}   "
              f"overlaps[0:8]={np.round(ov[:8],3)}",flush=True)
    results[a]=dstars[-1]
    print()
print("SUMMARY — stable core dim d*(>.97) at finest refinement, vs a:")
for a,d in results.items():
    print(f"   a={a} (log c={2*a:.1f}): d* = {d}")
xs=np.array(list(results.keys())); ys=np.array(list(results.values()))
if len(xs)>=2:
    sl=np.polyfit(2*xs,ys,1)[0]
    print(f"\n   d* vs log c slope = {sl:.2f}   (≈0 ⟹ BOUNDED core, Sub-goal H tractable; large + ⟹ growing/infinite)")
print("READ: a flat/bounded d* across a = a genuine FIXED finite core the reweighting must neutralize (tractable target).")
print("      A d* growing ~linearly in log c = the CCM ~1.8 log c growth is REAL structure, not grid noise (harder).")
