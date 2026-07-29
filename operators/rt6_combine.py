#!/usr/bin/env python3
"""§421 (ATTACK, the C-C combination frontier): C-C (2020) proved archimedean positivity via Toeplitz control
of prolates on the orthogonal complement of the cutoff. The OPEN step is combining the PRIME places. Bridge
idea: prime = P diag(ρ_prime) Pᵀ, ρ_prime(r)=Σ Λ(n)/√n cos(r log n) = sum of pure cosines = naturally Toeplitz.
Test whether, in the archimedean eigenbasis (the instrument's prolate-analog), the PRIME form restricted to the
archimedean-CRITICAL subspace (small archpole eigenvalues = C-C's orthogonal complement of the cutoff) is
LOW-RANK / bounded — the structure C-C's Toeplitz method could absorb. If prime|crit is low-rank & its norm <
archmargin there, the combination is controllable by the SAME mechanism; if it's full-rank/large, the prime
places genuinely need new input (the §332 second space / twin-prime)."""
import numpy as np, l2_fejer as LF   # use Fejér (the clean O(1) lens, §419)
BAND=LF.BAND
print(f"§421 prime form in the archimedean eigenbasis, restricted to the critical (small-archpole) subspace.")
print(f"Fejér instrument. band{BAND}. 'crit' = archpole eigenvalues below median (C-C orthogonal-complement analog).")
print(f"{'s':>5} {'n':>3} {'archmargin':>11} {'dim_crit':>8} {'rank(prime|crit)':>16} {'||prime|crit||':>14} {'ratio/arch':>10}")
for s in np.arange(5.5,2.6,-0.4):
    cs=LF.probes(BAND,s)
    if len(cs)<5: continue
    f=LF.build(cs); AP=f['archpole']; PR=f['prime']
    la,U=np.linalg.eigh(AP)           # archimedean eigenbasis, ascending
    archm=float(la[0])
    # critical subspace = eigenvectors with archpole eigenvalue below median (the "weak/cutoff-complement" block)
    k=max(2,len(cs)//2)
    Uc=U[:,:k]                        # k smallest-archpole directions
    Pc=Uc.T@PR@Uc                     # prime restricted to critical subspace
    ev=np.linalg.eigvalsh(Pc)
    rank=int(np.sum(np.abs(ev)>1e-9*max(1,np.abs(ev).max())))
    nrm=float(np.max(np.abs(ev)))
    print(f"{s:5.1f} {len(cs):3d} {archm:+11.4e} {k:8d} {rank:16d} {nrm:14.4e} {nrm/archm:10.3f}",flush=True)
print("""
reading: rank(prime|crit) << dim_crit AND ||prime|crit|| bounded => prime is LOW-RANK on the cutoff-complement,
  exactly the structure C-C's Toeplitz/prolate control could absorb => the combination step is plausibly a
  finite-rank perturbation of the proven archimedean positivity (a concrete crack). Full-rank & ||.||~archmargin
  => prime genuinely fills the critical subspace, needing the §332 second space / arithmetic (twin) input.
""")
