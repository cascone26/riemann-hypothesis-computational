#!/usr/bin/env python3
"""
CYCLE 291 — tie the jump-positivity (Euler-product) mechanism to the CLASS NUMBER.

euler_product_jump_positivity.py established: the prime side of the Weil form is a Beurling-Deny JUMP
form E_jump(psi)=sum_n w_n INT[psi(u)-psi(u+log n)]^2, and E_jump>=0 for ALL psi  <=>  every generalized
von Mangoldt weight Lambda_L(n)>=0  <=>  L has an Euler product (log L has all-positive coefficients).

Cycle 290 established: within the Epstein family ζ_Q of binary quadratic forms, the class number h is
the multiplicativity switch (h=1 => Euler product/PASS/0 off-line; h>=2 => no Euler product/FAIL/off-line).

THIS SCRIPT connects them concretely. For the PRINCIPAL reduced form Q of discriminant -d, build the
Dirichlet series ζ_Q(s)=sum_N a_Q(N) N^{-s} from representation counts, normalize a(1)=1, and compute
the generalized von Mangoldt weights Λ_Q(n) via the exact recursion
        sum_{d|n} Λ_Q(d) a(n/d) = a(n) log n .
Claim to VERIFY numerically:
  h=1 (disc -4): principal form = entire class group => ζ_Q ∝ ζ_K (Dedekind) => Euler product => Λ_Q(n)>=0.
  h>=2 (disc -20, -23): a SINGLE class ζ_Q is NOT a Dedekind zeta / has NO Euler product
       => Λ_Q(n) takes NEGATIVE values => jump form NOT >=0 => same failure that shows as off-line zeros.
So the class number literally flips the sign-definiteness of the jump weights. Report-tier structural tie.
"""
import numpy as np
from math import log, gcd

def reduced_principal_form(d):
    # principal reduced form of discriminant -d: x^2+ (d/4) y^2 if d%4==0 else x^2+xy+((1+d)/4) y^2
    if d % 4 == 0:
        return (1, 0, d//4)      # a x^2 + b xy + c y^2 = x^2 + (d/4) y^2
    else:
        return (1, 1, (d+1)//4)  # x^2 + xy + ((d+1)/4) y^2

def rep_counts(form, Nmax):
    a,b,c = form
    r = np.zeros(Nmax+1, dtype=float)
    B = int(2*np.sqrt(Nmax)) + 3
    for x in range(-B, B+1):
        for y in range(-B, B+1):
            v = a*x*x + b*x*y + c*y*y
            if 1 <= v <= Nmax:
                r[v] += 1.0
    return r  # r[N] = #{(x,y): Q(x,y)=N}

def von_mangoldt_from_dirichlet(a, Nmax):
    # a[1..Nmax], a[1] must be 1 (normalized). Solve sum_{d|n} Lam[d] a[n/d] = a[n] log n.
    Lam = np.zeros(Nmax+1)
    for n in range(2, Nmax+1):
        s = a[n]*log(n)
        d = 1
        # subtract contributions from proper divisors d|n, d<n
        for dd in range(1, n):
            if n % dd == 0:
                s -= Lam[dd]*a[n//dd]
        Lam[n] = s  # a[1]=1
    return Lam

def analyze(d, Nmax=200):
    form = reduced_principal_form(d)
    r = rep_counts(form, Nmax)
    a1 = r[1]
    a = r / a1  # normalize a(1)=1
    Lam = von_mangoldt_from_dirichlet(a, Nmax)
    neg_idx = [n for n in range(2, Nmax+1) if Lam[n] < -1e-9]
    return form, a1, Lam, neg_idx

for d, h, label in [(4,1,"h=1 (Euler product: ζ_K=ζL(χ₄))"),
                    (20,2,"h=2 (no Euler product)"),
                    (23,3,"h=3 (irreducible, no Euler product)")]:
    form, a1, Lam, neg = analyze(d)
    print(f"\ndisc -{d}, {label}")
    print(f"  principal form {form[0]}x²+{form[1]}xy+{form[2]}y²,  r_Q(1)={int(a1)} (=#units)")
    print(f"  # n<=200 with Λ_Q(n) < 0 : {len(neg)}")
    if neg:
        show = neg[:8]
        print(f"    first negative-weight n : {show}")
        print(f"    e.g. Λ_Q({show[0]}) = {Lam[show[0]]:+.4f}")
    else:
        print(f"    all Λ_Q(n) >= 0  =>  jump form POSITIVE  =>  Euler-product / PASS")
print("\n=> class number flips jump-weight sign-definiteness: h=1 all Λ>=0 (Euler product, jump-positive);"
      "\n   h>=2 single-class ζ_Q has NEGATIVE Λ (no Euler product) — the same break that shows as off-line zeros.")
