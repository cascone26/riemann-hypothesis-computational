#!/usr/bin/env python3
"""§413 (ATTACK): the instrument-INDEPENDENT fingerprint of the edge-closing.
§412 verify-catch: only the density-edge (margin->0 as probe count n -> RvM zero count N) is bandwidth-robust.
Now test: does λmin(B_analytic) close with a UNIVERSAL law across w when plotted vs the fill ratio n/N_RvM?
If λmin/λmin(deepest) collapses onto ONE curve across w=(4,3,2,1.5), the closing is instrument-independent
and its shape/exponent is a constant the second space must reproduce. If curves DON'T collapse, the closing
rate is also instrument-filtered and the only robust fact is the location (n=N) not the rate."""
import numpy as np, l2_edge_sweep as L
BAND=(14.0,44.0)
# RvM effective zero count in band (Riemann-von Mangoldt), independent of instrument:
def Nrvm(a,b):
    f=lambda T: (T/(2*np.pi))*np.log(T/(2*np.pi))-T/(2*np.pi)
    return f(b)-f(a)
N=Nrvm(*BAND); print(f"§413 edge-closing universality. band{BAND} RvM N={N:.2f} zeros")
print(f"{'w':>4} {'s':>5} {'n':>3} {'n/N':>5} {'λmin':>11} {'λ/λmax_w':>10}")
curves={}
for w in (4.0,3.0,2.0,1.5):
    pts=[]
    for s in np.arange(6.5,3.4,-0.25):
        cs=L.probes(BAND,s)
        if len(cs)<4: continue
        lm=float(np.linalg.eigvalsh(L.build_forms(cs,w)['full'])[0])
        pts.append((s,len(cs),lm))
    lmax=max(p[2] for p in pts)
    curves[w]=[(p[1]/N, p[2]/lmax, p) for p in pts]
    for ratio,norm,p in curves[w]:
        print(f"{w:4.1f} {p[0]:5.2f} {p[1]:3d} {ratio:5.2f} {p[2]:+11.3e} {norm:10.3e}",flush=True)
    print()
# collapse test: at matched n/N bins, is normalized λ the same across w?
print("COLLAPSE TEST — normalized λ/λmax at matched fill ratio n/N across w:")
print(f"{'n/N~':>6} " + "  ".join(f"w={w}" for w in (4.0,3.0,2.0,1.5)))
for target in (0.7,0.8,0.9,1.0,1.1):
    row=[]
    for w in (4.0,3.0,2.0,1.5):
        best=min(curves[w], key=lambda c: abs(c[0]-target))
        row.append(best[1] if abs(best[0]-target)<0.15 else np.nan)
    print(f"{target:6.2f} " + "  ".join(f"{r:8.2e}" if not np.isnan(r) else f"{'--':>8}" for r in row))
print("\nreading: rows ~constant across w => closing law is instrument-INDEPENDENT (universal edge fingerprint,")
print("its exponent/shape = a constant the 2nd space must match). rows vary wildly => only the LOCATION n=N is")
print("robust, the rate is filtered — then the invariant is purely the counting identity rank<=N (§404).")
