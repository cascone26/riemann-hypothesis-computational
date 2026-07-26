#!/usr/bin/env python3
"""§317b: corner law — J^{d,0} min/mean root gap vs d over a wide range (n=0 needs only d moments; cache=130).
Two-point §317 evidence gave exponent -1.5; nail it across d=16..120."""
import mpmath as mp, pickle, math
mp.mp.dps=50
gam=[mp.mpf(g) for g in pickle.load(open("operators/_gam130.pkl","rb"))]
def stat(d):
    co=[mp.binomial(d,k)*gam[k] for k in range(d,-1,-1)]
    scale=mp.sqrt(abs(co[0])*abs(co[-1])); co=[c/scale for c in co]
    rts=mp.polyroots(co,maxsteps=500,extraprec=260)
    im=max(abs(mp.im(r)) for r in rts)/max(abs(r) for r in rts)
    xs=sorted(mp.re(r) for r in rts)
    gaps=[xs[i+1]-xs[i] for i in range(len(xs)-1)]
    mean=(xs[-1]-xs[0])/(len(xs)-1)
    return min(gaps)/mean, im
prev=None
print("CORNER LAW: J^{d,0} min/mean gap vs d:",flush=True)
for d in [16,25,36,50,70,90,110,120]:
    g,im=stat(d)
    s=""
    if prev:
        pd,pg=prev
        s=f"  local exp={math.log(float(g/pg))/math.log(d/pd):+.3f}"
    print(f"  d={d:3d}: min/mean={mp.nstr(g,5)}  hyp={mp.nstr(im,3)}{s}",flush=True)
    prev=(d,g)
print("DONE",flush=True)
