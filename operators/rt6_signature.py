#!/usr/bin/env python3
"""§392 (angle 6): Hodge-index signature census — on LOW bands containing r=0 (never probed: the pole
has weight ONLY there, §358). Inertia of K=arch−prime and B_full=pole+K across fills on [0,B]. Question:
does the pole-anchored form show signature-(k+, rest) Hodge-index-type structure at/past criticality?"""
import numpy as np
import l2_edge_sweep as L
w = 3.0
def inertia(M, tol=1e-10):
    ev = np.linalg.eigvalsh(M); s = tol*max(1, abs(ev).max())
    return int((ev > s).sum()), int((ev < -s).sum()), int((abs(ev) <= s).sum())
print("LOW-BAND census (bands containing r=0; pole weighted; zeros start 14.13):")
print(f"{'band':>8} {'n':>3} {'inertia B_full':>14} {'inertia K=A-P':>14} {'lmin(B)':>10} {'pole_w':>9}")
for B in (20.0, 30.0, 40.0):
    for n in (4, 6, 8, 10, 12, 15):
        cs = np.linspace(0.5, B, n)
        f = L.build_forms(cs, w)
        z = 0.5j
        fp = np.exp(-((z-cs)**2)/(2*w*w)) + np.exp(-((z+cs)**2)/(2*w*w))
        pole = 2.0*np.real(np.outer(fp, fp))
        K = f['full'] - pole
        iB = inertia(f['full']); iK = inertia(K)
        lm = float(np.linalg.eigvalsh(f['full'])[0])
        pw = float(np.linalg.norm(pole, 2))
        print(f"[0.5,{B:4.0f}] {n:3d} {str(iB):>14} {str(iK):>14} {lm:+10.2e} {pw:9.2e}", flush=True)
print("\nreading: B_full PSD everywhere (RH true, on-line data). Hodge-index shape would be K with")
print("inertia (small+, many-) i.e. prime side dominating arch NEGATIVELY off a small positive part.")
print("If K stays mostly-positive even at r~0, the Hodge frame doesn't bite on this instrument; if K")
print("flips to (1-2)+/rest-negative at criticality on low bands, the one-sided certificate target opens.")
