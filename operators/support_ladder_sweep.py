#!/usr/bin/env python3
"""Empirical SUPPORT LADDER: Weil-form positivity margin vs test-function bandwidth.
Fejer probes: products F_jF_k have u-support |u| < 4W, so only primes with log n < 4W enter.
Rungs: 4W < log2=0.693 (archimedean only, proven positive - Connes-Consani rung), then p=2 enters
at 4W=log2, p=3 at log3=1.099, 4=2log2=1.386, 5=1.609, 7=1.946, 8,9, ...
Measure: min-eig of full = archpole - prime, normalized by ||archpole||, vs 4W.
Reading: where does the margin COLLAPSE from the (positive) archimedean rung toward the a=1 razor edge?
Which prime rung does the damage? This is the empirical map of the ladder our Track-2(a) would climb.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import l2_fejer as fj

BAND = (14.0, 44.0); SP = 1.5
print(f"{'4W':>6} {'primes in':>10} {'mineig(full)':>13} {'mineig/||AP||':>14} {'mineig(AP)':>12}")
for fourW in [0.5, 0.65, 0.72, 0.95, 1.05, 1.15, 1.30, 1.42, 1.55, 1.65, 1.90, 2.00, 2.30, 2.90, 3.50, 4.40]:
    fj.W = fourW/4.0
    cs = fj.probes(BAND, SP)
    f = fj.build(cs)
    FU = 0.5*(f['full']+f['full'].T); AP = 0.5*(f['archpole']+f['archpole'].T)
    me = float(np.linalg.eigvalsh(FU)[0]); meAP = float(np.linalg.eigvalsh(AP)[0])
    nAP = np.linalg.norm(AP)
    prim = [n for n in [2,3,4,5,7,8,9,11,13,16] if np.log(n) < fourW]
    print(f"{fourW:>6.2f} {str(prim)[:10]:>10} {me:>13.3e} {me/nAP:>14.3e} {meAP:>12.3e}", flush=True)
print("\nReading: margin at 4W<log2 = archimedean-only rung (should be ~PSD with margin, the proven rung);")
print("watch the normalized margin fall as each prime enters; the ladder's 'hard rung' = where it first")
print("hits ~0 (the razor edge reached with finitely many primes).")
