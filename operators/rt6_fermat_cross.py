#!/usr/bin/env python3
"""§425 (ATTACK, novel δ_p route, exploratory): the research (§424) flagged the Buium cross-prime δ_p∘δ_q
intersection invariant as GENUINELY UNEXPLORED. First correct probe: does the cross-prime Fermat-quotient
structure carry a definite / Hodge-type signature (the positivity the second space needs)? Fermat quotient
q_p(a) = ((a^{p-1}−1)/p) mod p ∈ {0..p−1} (Buium's δ_p at the base a, the standard bounded object). Build the
cross-prime correlation matrix G[p,q] over common bases a, normalized, and examine its signature. EXPLORATORY /
REPORT-tier: correct computation of a real object; the RH connection (does this seed the self-square pairing?)
is speculative and NOT claimed — just looking for suggestive structure vs a random-matrix null."""
import numpy as np
from sympy import primerange
def fermat_quotient_vec(p, bases):
    # q_p(a) = ((a^(p-1)-1)/p) mod p, for a coprime to p
    out=[]
    for a in bases:
        if a % p == 0: out.append(np.nan); continue
        val=(pow(a,p-1,p*p)-1)//p % p     # exact via mod p^2
        out.append(val/p)                 # normalize to [0,1)
    return np.array(out)
PRIMES=list(primerange(5,60))            # 5,7,...,59
bases=[a for a in range(2,200)]
V=[]
for p in PRIMES:
    v=fermat_quotient_vec(p,bases)
    V.append(v)
V=np.array(V)                            # rows = primes, cols = bases
# use only base columns valid for all primes (coprime to all) — drop any with nan
good=~np.any(np.isnan(V),axis=0); V=V[:,good]
# center each prime's vector (fluctuation about mean 1/2) and normalize
Vc=V-V.mean(axis=1,keepdims=True)
Vn=Vc/ (np.linalg.norm(Vc,axis=1,keepdims=True)+1e-12)
G=Vn@Vn.T                                # cross-prime correlation Gram
ev=np.linalg.eigvalsh(G)
def sig(e,tol=1e-9):
    m=max(1,abs(e).max()); return int((e>tol*m).sum()),int((abs(e)<=tol*m).sum()),int((e<-tol*m).sum())
print(f"§425 Fermat-quotient cross-prime correlation. primes={PRIMES[:8]}...(n={len(PRIMES)}), {V.shape[1]} bases")
print(f"eigenvalues (sorted): {np.array2string(ev,precision=3,max_line_width=200)}")
print(f"signature = {sig(G)}   trace={np.trace(G):.2f}  max|offdiag|={np.max(np.abs(G-np.diag(np.diag(G)))):.4f}")
# null comparison: random unit vectors same dims
rng=np.random.default_rng(0)
R=rng.standard_normal((len(PRIMES),V.shape[1])); R-=R.mean(1,keepdims=True); R/=np.linalg.norm(R,axis=1,keepdims=True)
Gr=R@R.T; evr=np.linalg.eigvalsh(Gr)
print(f"random-null eigenvalue spread: [{evr[0]:.3f}, {evr[-1]:.3f}]  vs Fermat: [{ev[0]:.3f}, {ev[-1]:.3f}]")
print(f"mean |offdiag| Fermat={np.mean(np.abs(G[~np.eye(len(PRIMES),dtype=bool)])):.4f} "
      f"vs random={np.mean(np.abs(Gr[~np.eye(len(PRIMES),dtype=bool)])):.4f}")
print("""
reading: if Fermat cross-prime correlations are INDISTINGUISHABLE from random (same eigenvalue spread, same
  offdiag magnitude), the Fermat quotients are cross-prime 'independent' — no special intersection structure at
  this (naive) level; the §332 pairing must live in the non-commutative δ_p∘δ_q composition, not the values.
  If Fermat shows STRUCTURE (definite signature, low-rank, large coherent offdiag), that's a real suggestive lead.
""")
