#!/usr/bin/env python3
"""§456 (combination structure): decompose the prime form's mass across the ARCHIMEDEAN spectrum. archpole=UΛU^T;
prime'=U^T prime U. Diagonal prime'[i][i] = prime's mass in archpole-eigenmode i (eigenvalue Λ[i]). Plot vs Λ[i]:
if prime concentrates where archpole is SMALL (weak) → arch−prime threatened there (the RH crux); if where archpole
is LARGE → arch dominates (harmless). This precisely characterizes whether the C-C semi-local combination can work."""
import numpy as np, l2_edge_sweep as L
BAND=(14.0,44.0); w=3.0
print("§456 prime mass vs archimedean eigenvalue. Gaussian, band[14,44].")
for s in (4.5,4.0,3.5):
    cs=L.probes(BAND,s); f=L.build_forms(cs,w)
    AP=f['archpole']; PR=f['prime']
    lam,U=np.linalg.eigh(AP)                 # archpole eigenvalues (ascending) + eigenvectors
    Pp=U.T@PR@U                               # prime in archpole eigenbasis
    pmass=np.diag(Pp)                          # prime mass per archpole mode
    print(f"\ns={s} n={len(cs)}: archpole eigenvalues (small→large) and prime mass there:")
    print(f"  {'archpole λ':>12} {'prime mass':>11} {'ratio pm/λ':>11}")
    for i in range(len(cs)):
        print(f"  {lam[i]:12.4e} {pmass[i]:11.4e} {pmass[i]/lam[i] if lam[i]>1e-15 else np.inf:11.3f}")
print("""
reading: if prime mass is LARGEST where archpole λ is SMALLEST (ratio pm/λ ≫1 at small λ), the primes attack
exactly the weak archimedean directions ⟹ arch−prime is threatened there = the RH crux is a few weak modes. If
prime mass is spread/large at big λ (harmless, arch dominates), positivity is easy. The ratio pm/λ at the smallest
archpole eigenvalue = the local 'stress'; RH ⟺ it stays ≤1 (arch dominates prime even in its weakest direction),
which §417/421 showed it does NOT trivially (coupled). This is the precise combination obstruction, mode by mode.
""")
