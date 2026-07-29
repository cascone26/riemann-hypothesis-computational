#!/usr/bin/env python3
"""§436 (cross-instrument, new): both RH-instruments singularize at the RH-content (§435). Test if the Nyman-
Beurling Gram's NEAR-NULL direction (smallest-eigenvalue eigvec of G, indexed by n=2..N) carries arithmetic
structure — does it concentrate on primes / correlate with μ(n)? Compare conceptually to the Weil critical
direction (§417: prime near-kernel). If both instruments' degenerate modes are prime/μ-structured, that's a
cross-instrument localization of where RH-content lives."""
import numpy as np
from sympy import mobius, isprime
J=400000; j=np.arange(1,J+1,dtype=np.float64); w=1.0/(j*(j+1))
N=60
V=np.array([(j%n)/n for n in range(2,N+1)]); Vw=V*w[None,:]; G=Vw@V.T
ev,evec=np.linalg.eigh(G)
ns=np.arange(2,N+1)
null=evec[:,0]                     # smallest-eigenvalue eigenvector (near-null direction)
mu=np.array([int(mobius(n)) for n in ns])
prime_mask=np.array([isprime(n) for n in ns])
# correlations
def corr(a,b): 
    a=a-a.mean(); b=b-b.mean(); return float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)+1e-12))
print(f"§436 Nyman-Beurling Gram near-null direction structure. N={N}, λmin(G)={ev[0]:.3e}, cond={ev[-1]/ev[0]:.1e}")
print(f"|null| weight on PRIMES vs composites: prime={np.sum(null[prime_mask]**2):.3f}  comp={np.sum(null[~prime_mask]**2):.3f}  (#prime={prime_mask.sum()},#comp={(~prime_mask).sum()})")
print(f"corr(null, μ(n)) = {corr(null,mu):+.3f}")
print(f"corr(null, 1/n) = {corr(null,1.0/ns):+.3f}   corr(null, Λ-like 1_prime) = {corr(null,prime_mask.astype(float)):+.3f}")
# top components by |null|
order=np.argsort(-np.abs(null))
print("top-10 |null| components (n, null_n, μ(n), prime?):")
for i in order[:10]:
    print(f"   n={ns[i]:3d}: {null[i]:+.3f}  μ={mu[i]:+d}  {'prime' if prime_mask[i] else 'comp'}")
print("""
reading: if the near-null direction concentrates on PRIMES / correlates with μ(n), the Nyman-Beurling degeneracy
carries the SAME arithmetic content as the Weil critical direction (§417 prime near-kernel) — cross-instrument
localization: RH-content lives in the prime/Möbius structure of both instruments' degenerate edge. If it's
featureless (spread, no μ/prime correlation), the two instruments' degeneracies are unrelated.
""")
