#!/usr/bin/env python3
"""§417B: the crux is ALIGNMENT (§417A: primes have 292x the energy to violate; only near-orthogonality of the
critical direction to the prime modes saves it). Characterize that orthogonality. v* = critical eigenvector
(bottom of archpole = where density margin vanishes). Prime form PR = Σ λ_i u_i u_i^T. Compute the overlaps
|⟨v*, u_i⟩| ranked by prime eigenvalue λ_i. Is v* orthogonal to the TOP prime modes (few-dim structure) and
only overlaps weak ones? Also cross-check w-robustness (is the near-orthogonality bandwidth-stable, unlike
§412's per-prime stiffness which was an artifact?)."""
import numpy as np, l2_edge_sweep as L
BAND=(14.0,44.0)
for w in (3.0,2.0,4.0):
    print(f"\n=== w={w} band{BAND} ===  |⟨v*, prime-eigvec_i⟩| by descending prime eigenvalue")
    for s in (4.5,4.0,3.6):
        cs=L.probes(BAND,s)
        if len(cs)<4: continue
        f=L.build_forms(cs,w); AP=f['archpole']; PR=f['prime']
        va, va_vec=np.linalg.eigh(AP); vstar=va_vec[:,0]        # bottom of archpole = critical dir
        pv,pu=np.linalg.eigh(PR)                                 # prime eigenpairs ascending
        order=np.argsort(pv)[::-1]                               # descending
        ov=[abs(float(vstar@pu[:,i])) for i in order]
        pvv=[pv[i] for i in order]
        # energy of v* explained by top-k prime modes
        cum=np.cumsum([o**2 for o in ov])
        print(f"  s={s} n={len(cs)}: primeλ_top={pvv[0]:.3f} overlaps=["+
              ", ".join(f"{o:.3f}" for o in ov[:6])+"]  top3-energy={:.3f}".format(cum[2] if len(cum)>2 else cum[-1]))
print("""
reading: if |⟨v*,u_top⟩| ~ 0 (v* avoids the biggest prime modes) and v*'s energy sits in the SMALL prime modes,
  the cancellation is structured: the density-critical direction lives where the primes are weak. If that holds
  across w (unlike the §412 artifact), it's a real geometric fact = RH-content: the archimedean-weak subspace
  and the prime-strong subspace are nearly orthogonal. The DIMENSION of the overlap = size of the real crux.
""")
