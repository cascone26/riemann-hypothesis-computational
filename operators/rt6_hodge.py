#!/usr/bin/env python3
"""§414 (ATTACK, toward the forced 2nd space): is B_analytic = pole+(arch−prime) a HODGE-INDEX inequality?
Weil's F_q proof: NS lattice has signature (1,ρ−1); the ample anchor (diagonal Δ) is the single + direction;
for classes ⟂ anchor, self-intersection ≤ 0 (Hodge index); Weil positivity = anchor dominates the neg primitive.
Translate to the explicit formula: candidate anchor = the rank-1 POLE term (ζ's s=1 pole, structurally +).
TEST: whiten by overlap G (Euclidean metric on probe space), split pole's positive direction (the anchor),
project (arch−prime) onto the anchor's orthogonal complement, and ask: is the projected primitive part
NEGATIVE-semidefinite? If yes => the explicit formula IS a Hodge-index inequality (RH ⟺ anchor dominates);
if indefinite => there are + primitive directions the anchor must ALSO dominate (not a clean Hodge structure)."""
import numpy as np, l2_edge_sweep as L
BAND=(14.0,44.0)
def whiten(M,G,rc=1e-12):
    s,V=np.linalg.eigh(G); keep=s>rc*s[-1]; Vs=V[:,keep]/np.sqrt(s[keep])
    return Vs.T@M@Vs
print(f"§414 Hodge-index recast. band{BAND} w=3. anchor=pole (rank-1 +). is (arch−prime)⟂anchor neg-semidef?")
print(f"{'s':>5} {'n':>3} sig(pole)   sig(arch−prime)     anchor_λ   prim⟂ sig(+/0/−)   prim⟂ λmax  full≥0?")
for s in np.arange(6.5,3.4,-0.5):
    cs=L.probes(BAND,s)
    if len(cs)<4: continue
    f=L.build_forms(cs,3.0); G=f['overlap']
    pole=f['archpole']-(f['archpole']-2.0*np.real(np.outer(  # recompute pole cleanly
        np.exp(-((0.5j-cs)**2)/(18.0))+np.exp(-((0.5j+cs)**2)/(18.0)),
        np.exp(-((0.5j-cs)**2)/(18.0))+np.exp(-((0.5j+cs)**2)/(18.0)))))
    # cleaner: recompute pole directly
    fp=np.exp(-((0.5j-cs)**2)/(2*9.0))+np.exp(-((0.5j+cs)**2)/(2*9.0))
    pole=2.0*np.real(np.outer(fp,fp))
    arch=f['archpole']-pole
    R=arch-f['prime']           # the primitive candidate (arch − prime)
    Pw=whiten(pole,G); Rw=whiten(R,G); Fw=whiten(f['full'],G)
    # anchor = top eigvec of whitened pole
    ev,evec=np.linalg.eigh(Pw); anchor=evec[:,-1]
    Q=np.eye(len(anchor))-np.outer(anchor,anchor)   # projector ⟂ anchor
    Rperp=Q@Rw@Q
    def sig(M,tol=1e-9):
        e=np.linalg.eigvalsh(M); m=max(1.0,abs(e).max())
        return (int(np.sum(e>tol*m)),int(np.sum(np.abs(e)<=tol*m)),int(np.sum(e<-tol*m)))
    sp=sig(Pw); sr=sig(Rw); sperp=sig(Rperp)
    rperp_ev=np.linalg.eigvalsh(Rperp)
    lmax_perp=float(rperp_ev[-1]) if len(rperp_ev) else 0.0
    fmin=float(np.linalg.eigvalsh(Fw)[0])
    print(f"{s:5.1f} {len(cs):3d}  {str(sp):>9} {str(sr):>12}   {ev[-1]:+.2e}   {str(sperp):>12}   "
          f"{lmax_perp:+.2e}  {'YES' if fmin>=-1e-9 else 'no'}",flush=True)
print("""
reading: sig=(#+,#0,#−). If prim⟂ sig has ZERO '+' (all 0/−), (arch−prime) is neg-semidef ⟂ anchor =>
  Weil positivity IS Hodge-index (RH ⟺ single ample anchor dominates a purely negative primitive part) —
  a genuine structural bridge to the F_q proof. If prim⟂ still has '+' directions (lmax_perp>0), the anchor
  is NOT the whole story: those + primitive modes mean the form is NOT a one-anchor Hodge inequality here.
""")
