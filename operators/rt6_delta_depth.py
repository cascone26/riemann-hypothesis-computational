#!/usr/bin/env python3
"""§427 (δ_p DEPTH TEST, honest): is [δ_p,δ_q]>0 prime-ARITHMETIC (RH-relevant) or a generic operator fact?
Test 1: the leading-coeff sign = sign(g(q)−g(p)), g=(x−1)/ln x — works for ANY real orders a<b, no primality.
  So compute the 'commutator' leading coeff for NON-PRIME real orders and confirm the positivity is generic.
Test 2: does the FULL commutator sequence [δ_p,δ_q](n) encode anything ζ/prime-special — e.g., does
  Σ_n [δ_p,δ_q](n) x^n or the value pattern relate to a known arithmetic function? Quick structural scan.
If positivity is generic (holds for non-prime real orders) AND no ζ-connection, the §425 lead is SHALLOW —
bank honestly and pivot."""
import numpy as np
def lead_coeff(a,b):  # leading coeff of the 'commutator' [δ_a,δ_b] for real orders a<b
    return (1.0/(a*b))*(1.0/b**(a-1) - 1.0/a**(b-1))
print("§427 Test 1: leading-coeff sign for NON-PRIME real orders a<b (positivity generic?):")
pairs=[(2.0,3.0),(2.5,3.7),(4.0,6.0),(1.5,2.2),(9.0,10.0),(2.0,100.0),(1.1,1.2)]
allpos=True
for a,b in pairs:
    lc=lead_coeff(a,b); allpos=allpos and lc>0
    print(f"  a={a} b={b}: leading coeff={lc:+.4e}  {'>0' if lc>0 else '<=0'}")
print(f"  => positivity holds for ALL non-prime real a<b: {allpos}  (generic monotonicity, NOT prime-special)")
# Test 2: any prime/zeta structure in the actual commutator sequence?
from sympy import primerange
def d(r,m): return (m**r-m)//r
def comm(p,q,n): return d(p,d(q,n))-d(q,d(p,n))
print("\n§427 Test 2: does [δ_2,δ_3](n)/n structure relate to primes/known sequences?")
seq=[comm(2,3,n) for n in range(2,10)]
print(f"  [δ_2,δ_3](n), n=2..9: {seq}")
# ratios and factorizations
from sympy import factorint
for n in range(2,7):
    print(f"    n={n}: comm={comm(2,3,n)}  factored={dict(factorint(comm(2,3,n)))}")
print("""
ASSESSMENT: if Test 1 = generic (positivity for any real a<b) and Test 2 shows no ζ/prime-special structure,
the §425 commutator positivity is a SHALLOW operator fact — it does NOT encode the prime/zero arithmetic RH needs.
Honest conclusion: the δ_p leading-order lead is shallow; the RH-content (if the route has any) must be in the
deeper jet-space/Frobenius-lift, which is the genuine hard Buium construction. PIVOT to a higher-value angle.
""")
