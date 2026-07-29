#!/usr/bin/env python3
"""§414b: the §414 anchor (pole) was numerically ZERO in band[14,44] (far from s=1). The ζ-pole term
F(i/2) is O(1) only for probes near γ~0. Redo the Hodge-index test on a WIDE band (2..44) that INCLUDES
the pole region (low probes see the anchor) AND zeros (high probes). Question: with the anchor VISIBLE,
does (arch−prime) become neg-semidef ⟂ anchor (=> Hodge-index structure), or does full stay 'positive
because arch−prime is positive' (=> no anchor structure, just Weil positivity)?"""
import numpy as np, l2_edge_sweep as L
def whiten(M,G,rc=1e-12):
    s,V=np.linalg.eigh(G); keep=s>rc*s[-1]; Vs=V[:,keep]/np.sqrt(s[keep])
    return Vs.T@M@Vs
def poleterm(cs,w):
    fp=np.exp(-((0.5j-cs)**2)/(2*w*w))+np.exp(-((0.5j+cs)**2)/(2*w*w)); return 2.0*np.real(np.outer(fp,fp))
def sig(M,tol=1e-9):
    e=np.linalg.eigvalsh(M); m=max(1.0,abs(e).max())
    return (int(np.sum(e>tol*m)),int(np.sum(np.abs(e)<=tol*m)),int(np.sum(e<-tol*m)))
w=3.0
# first: how big is the pole term diag vs arch/prime, across probe center?
print("pole-term magnitude vs probe center c (w=3):  ||pole col|| relative to ||arch−prime col||")
for c in (2.0,5.0,8.0,14.0,25.0,40.0):
    cs=np.array([c]); fp=np.exp(-((0.5j-cs)**2)/(2*w*w))+np.exp(-((0.5j+cs)**2)/(2*w*w))
    print(f"  c={c:5.1f}: pole_diag={2*np.real(fp*fp)[0]:.3e}")
print()
BAND=(2.0,44.0)
print(f"§414b Hodge test on WIDE band{BAND} w={w} (anchor now visible at low probes).")
print(f"{'s':>5} {'n':>3}  anchor_λ    sig(arch−prime)   prim⟂ sig(+/0/−)  prim⟂ λmax  full λmin  full≥0?")
for s in np.arange(6.5,2.9,-0.5):
    cs=L.probes(BAND,s)
    if len(cs)<4: continue
    f=L.build_forms(cs,w); G=f['overlap']
    pole=poleterm(cs,w); arch=f['archpole']-pole; R=arch-f['prime']
    Pw=whiten(pole,G); Rw=whiten(R,G); Fw=whiten(f['full'],G)
    ev,evec=np.linalg.eigh(Pw); anchor=evec[:,-1]
    Q=np.eye(len(anchor))-np.outer(anchor,anchor); Rperp=Q@Rw@Q
    rp=np.linalg.eigvalsh(Rperp)
    print(f"{s:5.1f} {len(cs):3d}  {ev[-1]:+.2e}  {str(sig(Rw)):>12}   {str(sig(Rperp)):>12}   "
          f"{rp[-1]:+.2e}  {float(np.linalg.eigvalsh(Fw)[0]):+.2e}  "
          f"{'YES' if np.linalg.eigvalsh(Fw)[0]>=-1e-9 else 'no'}",flush=True)
print("\nreading: anchor_λ now O(1)? then the test is meaningful. prim⟂ ZERO '+' => Hodge-index structure")
print("appears once the anchor is in view. prim⟂ still '+' => no single-anchor Hodge form even with pole visible.")
