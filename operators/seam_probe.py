#!/usr/bin/env python3
"""§317: SEAM PROBE. The uniformity crux is the matched-asymptotics boundary n ~ N_dir(d)=0.745 d^{3/2}
between the Hermite (outer) and Xi (inner) regimes. Measure the hyperbolicity margin PROFILE across the
seam: normalized min root gap of J^{d,n} (min gap / mean gap) as n runs 0..105, at d=16 (seam ~48) and
d=25 (seam ~93) — both seams inside our 130-moment reach. Reference: Hermite He_d's own min/mean gap
(the n->inf limit). Question: does the margin DIP at/inside the seam (tension building toward the Xi
corner) or improve monotonically with n? Charisma: HANDLE (measured profile)."""
import mpmath as mp, pickle
mp.mp.dps=60
gam=[mp.mpf(g) for g in pickle.load(open("operators/_gam130.pkl","rb"))]
def roots_J(d,n):
    co=[mp.binomial(d,k)*gam[n+k] for k in range(d,-1,-1)]
    scale=mp.sqrt(abs(co[0])*abs(co[-1])); co=[c/scale for c in co]
    rts=mp.polyroots(co,maxsteps=400,extraprec=250)
    xs=sorted(mp.re(r) for r in rts)
    return xs
def gapstat(xs):
    d=len(xs)-1
    gaps=[xs[i+1]-xs[i] for i in range(len(xs)-1)]
    mean=(xs[-1]-xs[0])/ (len(xs)-1)
    return min(gaps)/mean
def hermite_ref(d):
    P=[mp.mpf(1)]; Q=[mp.mpf(0),mp.mpf(1)]
    for k in range(1,d):
        R=[mp.mpf(0)]+Q
        for i in range(len(P)): R[i]-=k*P[i]
        P,Q=Q,R
    rts=mp.polyroots(list(reversed(Q)),maxsteps=300,extraprec=150)
    return gapstat(sorted(mp.re(r) for r in rts))
for d in (16,25):
    ref=hermite_ref(d)
    seam=int(0.745*d**1.5)
    print(f"d={d}: Hermite ref min/mean gap = {mp.nstr(ref,5)}; seam N_dir ~ {seam}",flush=True)
    print("   n : min/mean gap   (vs ref)",flush=True)
    for n in list(range(0,106,5)):
        if n+d>130: break
        g=gapstat(roots_J(d,n))
        mark="<-- seam" if abs(n-seam)<=2 else ""
        print(f"  {n:3d}: {mp.nstr(g,5)}   ({mp.nstr(g/ref,4)}x ref) {mark}",flush=True)
    print("",flush=True)
print("DONE",flush=True)
