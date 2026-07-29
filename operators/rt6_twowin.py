#!/usr/bin/env python3
"""§393 (angle 3): two-window entanglement — joint Gram on disjoint bands vs direct sum.
Cross-block norm, joint lmin vs min(individual), joint numerical rank vs sum, vs separation+fill."""
import numpy as np
from scipy.special import digamma
import l2_edge_sweep as L
w = 3.0
L.RGRID = np.arange(-250.0, 250.0, 0.002)
L.OMEGA = np.real(digamma(0.25 + 0.5j*L.RGRID)) - np.log(np.pi)
def nrank(M, tol=1e-9):
    ev = np.linalg.eigvalsh(M); return int((ev > tol*max(1e-300, ev[-1])).sum())
print(f"{'W2start':>7} {'fill':>5} {'lmin1':>9} {'lmin2':>9} {'lminJ':>9} {'J/min':>7} {'rk1+rk2':>7} {'rkJ':>4} {'||X||':>9}")
for W2s in (50.0, 100.0, 170.0):
    for s in (5.0, 4.0, 3.5):
        cs1 = np.arange(14.0, 44.0+1e-9, s); cs2 = np.arange(W2s, W2s+30.0+1e-9, s)
        csJ = np.concatenate([cs1, cs2])
        f1 = L.build_forms(cs1, w); f2 = L.build_forms(cs2, w); fJ = L.build_forms(csJ, w)
        l1 = float(np.linalg.eigvalsh(f1['full'])[0]); l2 = float(np.linalg.eigvalsh(f2['full'])[0])
        lJ = float(np.linalg.eigvalsh(fJ['full'])[0])
        n1 = len(cs1)
        X = fJ['full'][:n1, n1:]
        rk = nrank(fJ['full']); r12 = nrank(f1['full']) + nrank(f2['full'])
        print(f"{W2s:7.0f} {s:5.1f} {l1:9.2e} {l2:9.2e} {lJ:9.2e} {lJ/min(l1,l2):7.3f} {r12:7d} {rk:4d} "
              f"{np.linalg.norm(X,2):9.2e}", flush=True)
print("\nreading: J/min<<1 or rkJ<rk1+rk2 at subcritical fill = cross-window entanglement (global")
print("interpolation rigidity); J/min~1 and additive ranks = windows independent (certificates local).")
