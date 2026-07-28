#!/usr/bin/env python3
"""§388a: reconstruct the RT's connected-6-point local cumulant by matching their checkpoint:
p=5, s=36, five shape values {+0.60, +3.49, +3.49, -7.56, -7.56}.
Candidate object: local ALTERNATING (connected/M-S S0-style) cumulant at prime p for a
configuration X of distinct positions:  A_p(X) = sum_{J subseteq X} (-1)^{|X|-|J|} w_p(J),
w_p(J) = (1 - nu_p(J)/p)/(1-1/p)^{|J|},  nu_p = #distinct residues mod p  (w_p(empty)=1).
Parametrizations scanned (all at s=36):
  P1: 4-point gap-triple  X = {0, d1, d1+d2, 36},  (d1,d2,d3) even, sum 36.
  P2: 6-point two-triple  X = {0,d1,d1+d2} u {36+0, 36+e1, 36+e1+e2}? (shifted copy, d=e)
  P3: 6-point palindrome  X = {0, d1, d1+d2, 36-d1-d2, 36-d1, 36} (triple + reflected triple)
Print A_5 for all shapes in each parametrization; look for the five-value pattern.
"""
import numpy as np
from itertools import combinations

def w_p(J, p):
    if not J: return 1.0
    nu = len(set(x % p for x in J))
    return (1 - nu / p) / (1 - 1 / p) ** len(J)

def A_p(X, p):
    X = sorted(set(X))
    tot = 0.0
    n = len(X)
    for k in range(n + 1):
        for J in combinations(X, k):
            tot += (-1) ** (n - k) * w_p(J, p)
    return tot

s = 36
p = 5
target = sorted([0.60, 3.49, 3.49, -7.56, -7.56])

def gap_triples(s):
    out = []
    for d1 in range(2, s - 3, 2):
        for d2 in range(2, s - d1 - 1, 2):
            d3 = s - d1 - d2
            if d3 >= 2:
                out.append((d1, d2, d3))
    # unordered shapes
    return sorted(set(tuple(sorted(t)) for t in out))

shapes = gap_triples(s)
print(f"P1: 4-point gap-triple {{0,d1,d1+d2,36}}, p=5 — {len(shapes)} unordered shapes")
vals1 = {}
for sh in shapes:
    d1, d2, d3 = sh
    v = A_p([0, d1, d1 + d2, s], p)
    vals1[sh] = v
# distinct values and counts
from collections import Counter
c1 = Counter(round(v, 2) for v in vals1.values())
print("  value histogram:", dict(sorted(c1.items())))
hits = [sh for sh, v in vals1.items() if abs(v - 0.60) < 0.02 or abs(v - 3.49) < 0.02 or abs(v + 7.56) < 0.02]
print("  shapes near target values:", [(sh, round(vals1[sh], 3)) for sh in hits][:12])

print(f"\nP3: 6-point palindrome {{0,d1,d1+d2,36-d1-d2,36-d1,36}}, p=5")
vals3 = {}
for sh in shapes:
    for perm in set([(sh[0], sh[1], sh[2]), (sh[0], sh[2], sh[1]), (sh[1], sh[0], sh[2])]):
        d1, d2, d3 = perm
        X = [0, d1, d1 + d2, s - d1 - d2, s - d1, s]
        if len(set(X)) == 6:
            vals3[perm] = A_p(X, p)
c3 = Counter(round(v, 2) for v in vals3.values())
print("  value histogram:", dict(sorted(c3.items())))
hits3 = [sh for sh, v in vals3.items() if abs(v - 0.60) < 0.03 or abs(v - 3.49) < 0.03 or abs(v + 7.56) < 0.03]
print("  shapes near target:", [(sh, round(vals3[sh], 3)) for sh in hits3][:12])

# P4: maybe the five values come from the five RESIDUE CLASSES of a single parameter mod 5.
print("\nP4: A_5 as function of d mod 5 for symmetric 3-term progression {0,d,2d,...}? quick scans:")
for d in range(2, 12, 2):
    X4 = [0, d, 36 - d, 36]
    print(f"  X={{0,{d},{36-d},36}}: A_5 = {A_p(X4, p):+.3f}", end=";  ")
print()
print("\nreading: if any parametrization's value multiset contains the target 5-pattern")
print(f"  {target} on its tested cells, that is the RT's object — proceed to O-kernel build.")
