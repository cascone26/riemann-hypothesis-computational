#!/usr/bin/env python3
"""§415c (honest reframe): §415b showed off-diagonal is ARCH-dominated (prime only ~14%). So the real
structure of full=archpole−prime: archpole is the dominant form whose OWN λmin vanishes at the edge (overlap
degeneracy = RvM saturation); prime is a SMALL correction. RH⟺the small prime correction never pushes λmin<0.
Measure: λmin(archpole) vs λmin(full), and the prime 'push' = λmin(archpole)−λmin(full). Is the push SMALLER
than archpole's margin everywhere (RH holds with room) or does it become comparable (tight) at the edge?"""
import numpy as np, l2_edge_sweep as L
BAND=(14.0,44.0)
print(f"§415c archpole margin vs prime push. band{BAND} w=3.")
print(f"{'s':>5} {'n':>3} {'λmin(archpole)':>14} {'λmin(full)':>11} {'prime push':>11} {'push/archmargin':>15}")
for s in np.arange(6.5,2.9,-0.25):
    cs=L.probes(BAND,s)
    if len(cs)<4: continue
    f=L.build_forms(cs,3.0)
    la=float(np.linalg.eigvalsh(f['archpole'])[0])
    lf=float(np.linalg.eigvalsh(f['full'])[0])
    push=la-lf
    print(f"{s:5.1f} {len(cs):3d} {la:+14.4e} {lf:+11.4e} {push:+11.4e} {push/la if la>1e-15 else np.nan:15.3f}",flush=True)
print("""
reading: if λmin(archpole) itself → 0 at the edge, the density degeneracy (NOT primes) sets the scale, and
  the prime push is a small rider. push/archmargin < 1 everywhere => RH holds with room on this instrument
  (the primes never big enough to cross). push/archmargin → 1 at edge => tight (primes exactly fill the
  vanishing archimedean margin) — THAT tightness is the real RH-content, and it's between prime size and the
  archimedean DENSITY margin, not a pair-correlation. This is the honest 'what must the 2nd space control'.
""")
