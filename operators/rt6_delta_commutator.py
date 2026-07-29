#!/usr/bin/env python3
"""§425b (novel δ_p route, exploratory): the intersection content (if any) is in the NON-COMMUTATIVE δ_p∘δ_q
composition (§425: values are random). Compute Buium's exact commutator [δ_p,δ_q](n)=δ_pδ_q(n)−δ_qδ_p(n),
δ_r(m)=(m^r−m)/r (exact integers). Look for structure: is it 0? definite sign? a clean bilinear/divisibility
pattern in (p,q,n)? Buium's theory says these encode 'arithmetic curvature' — a nonzero, structured commutator
is the candidate seed for the missing self-square pairing. EXPLORATORY: correct exact computation, no RH claim."""
from sympy import primerange
def d(r,m): return (m**r - m)//r
def comm(p,q,n): return d(p,d(q,n)) - d(q,d(p,n))
PRIMES=list(primerange(2,12))   # 2,3,5,7,11
print("§425b Buium commutator [δ_p,δ_q](n) exact. Looking for sign/divisibility/bilinear structure.")
for (p,q) in [(2,3),(2,5),(3,5),(2,7),(3,7),(5,7)]:
    row=[]
    for n in range(2,9):
        c=comm(p,q,n)
        row.append(c)
    # structure diagnostics: signs, gcd, divisibility by pq
    signs=[('+' if x>0 else '-' if x<0 else '0') for x in row]
    from math import gcd
    g=0
    for x in row: g=gcd(g,abs(x))
    print(f" p={p} q={q}: signs={''.join(signs)}  gcd|comm|={g}  comm(n=2..4)={row[:3]}")
# antisymmetry check + a cross-prime pairing candidate: normalized comm growth exponent
print("\ngrowth: log|comm(p,q,n)| / log(n) should ~ pq (leading term n^{pq}); check the leading behavior:")
import numpy as np
for (p,q) in [(2,3),(3,5)]:
    ns=list(range(3,8)); lv=[np.log(abs(comm(p,q,n))) for n in ns]
    slope=np.polyfit([np.log(n) for n in ns],lv,1)[0]
    print(f" p={p} q={q}: empirical exponent={slope:.2f}  (pq={p*q})")
print("""
reading: definite sign + clean gcd (e.g. always divisible by pq or by a binomial factor) => the commutator is a
  STRUCTURED arithmetic invariant (Buium curvature) = a real candidate for the self-square pairing seed, worth
  a careful build. sign-random / gcd=1 / no pattern => even the non-commutative structure is generic at this
  level; the pairing needs the full jet-space/Frobenius-lift machinery, not the raw integer commutator.
""")
