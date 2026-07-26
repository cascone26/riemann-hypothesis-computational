#!/usr/bin/env python3
"""§339b: (k,R) frontier through the plunge floor, band [14,44], w=3.
The §338 chain (k=3, R=3) fails below s=3.9 due to ~5e-6 METHOD slack. What (k,R) certifies
s = 3.8 .. 3.5? Measures how the required structure grows as the margin plunges."""
import numpy as np
import l2_edge_sweep as L
from resonant_cert import prime_terms
from band_sweep_cert import gap_bound  # noqa  (same gap_bound signature)

w = 3.0
for s in (3.9, 3.8, 3.7, 3.6, 3.5):
    cs = L.probes(L.BAND, s)
    f = L.build_forms(cs, w)
    A = f['archpole']; P = f['prime']
    true_m = float(np.linalg.eigvalsh(f['full'])[0])
    terms = prime_terms(cs)
    order = sorted(terms, key=lambda n: -np.linalg.norm(terms[n], 2))
    nmax = len(cs)
    print(f"\ns={s}  n={nmax}  true margin {true_m:+.3e}   bound table (rows k, cols R):")
    Rlist = [3, 5, 8, 12]
    klist = [3, 4, 5, 6, min(7, nmax - 1)]
    print("      " + " ".join(f"R={R:<9d}" for R in Rlist))
    found = None
    for k in klist:
        vals = []
        for R in Rlist:
            Pres = sum(terms[n] for n in order[:R])
            b = gap_bound(A - Pres, P - Pres, k=k)
            vals.append(b)
            if b > 0 and found is None:
                found = (k, R)
        print(f"  k={k}: " + " ".join(f"{v:+.2e}" for v in vals))
    print(f"  first certifying (k,R): {found if found else 'NONE in grid'}")
