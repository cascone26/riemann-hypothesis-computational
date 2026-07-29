#!/usr/bin/env python3
"""§426 (develop §425 lead, exploratory): the sign-definite antisymmetric leading-coefficient cross-prime form
c(p,q)=(1/pq)[1/q^{p-1}−1/p^{q-1}] (p<q; antisymmetric). Build A[p,q]=c(p,q), examine its SYMPLECTIC spectrum
(eigenvalues ±iμ_k of a real antisymmetric matrix). Questions: (1) is A full-rank (nondegenerate symplectic form
= a genuine Kähler structure) or degenerate? (2) do the symplectic eigenvalues μ_k show structure (decay, prime
relation)? (3) Pfaffian sign. EXPLORATORY, no RH claim."""
import numpy as np
from sympy import primerange
PR=list(primerange(2,50))
n=len(PR)
def c(p,q):  # leading coeff of [δ_p,δ_q](·)/·^{pq}, defined for p<q
    return (1.0/(p*q))*(1.0/q**(p-1) - 1.0/p**(q-1))
A=np.zeros((n,n))
for i in range(n):
    for j in range(i+1,n):
        v=c(PR[i],PR[j]); A[i,j]=v; A[j,i]=-v
# symplectic spectrum: eigenvalues of A are ±i μ_k
ev=np.linalg.eigvals(A)
mu=np.sort(np.abs(ev.imag))[::-1]
mu=mu[mu>1e-14][:n//2+1]                      # positive symplectic eigenvalues
rank=np.linalg.matrix_rank(A,tol=1e-12)
print(f"§426 symplectic spectrum of the leading-coeff cross-prime form. {n} primes {PR[0]}..{PR[-1]}")
print(f"rank(A)={rank}/{n}  (full/even => nondegenerate symplectic = genuine Kähler candidate)")
print(f"top symplectic eigenvalues μ_k: {np.array2string(mu[:8],precision=4)}")
print(f"μ_k decay ratios: {np.array2string(mu[:6]/mu[1:7],precision=3)}")
# Pfaffian sign (via eigenvalues): sign = sign(prod of +imag parts)
posimag=np.sort(ev.imag[ev.imag>1e-14])
print(f"# nonzero symplectic pairs: {len(posimag)}  (n/2={n//2})")
# does log μ_k relate to primes? crude check vs 1/p^2 scale (dominant term ~ 1/(pq p^{q-1}) tiny for big q)
print(f"\nμ_k vs A's dominant entries |c(2,q)| (the largest, since 1/q^1 slow-decay):")
for q in PR[1:6]:
    print(f"  c(2,{q})={c(2,q):+.5f}", end="  ")
print()
print("""
reading: full-rank A => nondegenerate arithmetic symplectic/Kähler form on the prime lattice — a real structure.
  Structured μ_k (clean decay / prime-related) => worth pairing with a symmetric metric g to test a Hodge-index
  inequality (positivity of g relative to the Kähler form). Degenerate / featureless => the leading-coeff form is
  too crude; need the full δ_p-jet space. Either way this is the concrete next handle on the novel second-space route.
""")
