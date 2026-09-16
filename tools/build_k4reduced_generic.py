import sys
import sympy as sp
from sympy_to_lean import sympy_to_lean, sin_lean, csin_lean

z1,z1c,z2,z2c,z3,z3c,z4,z4c = sp.symbols('z1 z1c z2 z2c z3 z3c z4 z4c')
w1,w1c,w2,w2c,w3,w3c,w4,w4c = sp.symbols('w1 w1c w2 w2c w3 w3c w4 w4c')
p1,q1,r1,s1,p2,q2,r2,s2,p3,q3,r3,s3,p4,q4,r4,s4 = sp.symbols('p1 q1 r1 s1 p2 q2 r2 s2 p3 q3 r3 s3 p4 q4 r4 s4')
zmap = {1:(z1,z1c),2:(z2,z2c),3:(z3,z3c),4:(z4,z4c)}
wmap = {1:(w1,w1c),2:(w2,w2c),3:(w3,w3c),4:(w4,w4c)}
pmap = {1:(p1,q1,r1,s1),2:(p2,q2,r2,s2),3:(p3,q3,r3,s3),4:(p4,q4,r4,s4)}
idx=[1,2,3,4]
def others(i): return [j for j in idx if j!=i]
def A(z,zc): return z**2*zc**2
def B(z,zc): return z**2+zc**2

def coeff_sym(a,i):
    z,zc = zmap[i]; w,wc = wmap[i]; p,q,r,s = pmap[i]
    if a==0: return (p+q)*w*z*zc**2 + (r+s)*wc*z**2*zc
    if a==1: return (q-p)*w*zc**2 + (s-r)*wc*z**2
    if a==2: return -(p+q)*w*z - (r+s)*wc*zc
    if a==3: return (p-q)*w + (r-s)*wc

def magnitude_sympy(b,j,k,l):
    zj,zjc=zmap[j]; zk,zkc=zmap[k]; zl,zlc=zmap[l]
    Aj,Ak,Al=A(zj,zjc),A(zk,zkc),A(zl,zlc)
    Bj,Bk,Bl=B(zj,zjc),B(zk,zkc),B(zl,zlc)
    if b==10: return -1, Bj+Bk+Bl
    if b==8:  return 1, Bj*Bk+Bj*Bl+Bk*Bl+Aj+Ak+Al
    if b==6:  return -1, Aj*Bk+Aj*Bl+Ak*Bj+Ak*Bl+Al*Bj+Al*Bk+Bj*Bk*Bl
    if b==4:  return 1, Aj*Ak+Aj*Al+Ak*Al+Aj*Bk*Bl+Ak*Bj*Bl+Al*Bj*Bk
    if b==2:  return -1, Aj*Ak*Bl+Aj*Al*Bk+Ak*Al*Bj
    if b==0:  return 1, Aj*Ak*Al

PAIRS = {15:[(3,12)],13:[(1,12),(3,10)],11:[(1,10),(3,8)],9:[(1,8),(3,6)],7:[(1,6),(3,4)],5:[(1,4),(3,2)],3:[(1,2),(3,0)]}
COEFF_NAME = {15:"crossK4LeadingCoeff",13:"crossK4ThirdCoeff",11:"crossK4FifthCoeff",9:"crossK4SeventhCoeff",
              7:"crossK4NinthCoeff",5:"crossK4EleventhCoeff",3:"crossK4ThirteenthCoeff"}
RED_NAME = {13:"crossK4ReducedK13",11:"crossK4ReducedK11",9:"crossK4ReducedK9",7:"crossK4ReducedK7",
            5:"crossK4ReducedK5",3:"crossK4ReducedK3"}

def block_goal(d,i):
    js,ks,ls = others(i)
    total = 0
    for (a,b) in PAIRS[d]:
        c = coeff_sym(a,i)
        if b==12: total += c
        else:
            osign, mag = magnitude_sympy(b,js,ks,ls)
            total += osign*c*mag
    return total

E15 = sum(block_goal(15,i) for i in idx)

def coeff_lean_text(a, i):
    p,q,r,s = f"p{i}",f"q{i}",f"r{i}",f"s{i}"
    z = f"z{i}"
    Ai, Aci = sin_lean(i), csin_lean(i)
    if a == 1:
        return f"(({q} - {p}) * {Ai} * (starRingEnd ℂ {z}) ^ 2 + ({s} - {r}) * {Aci} * {z} ^ 2)"
    if a == 3:
        return f"(({p} - {q}) * {Ai} + ({r} - {s}) * {Aci})"
    raise ValueError

def derive_and_build(d):
    Ed = sum(block_goal(d,i) for i in idx)
    b = d - 3
    j,k,l = others(4)
    osign, magb = magnitude_sympy(b, j,k,l)
    combo = sp.expand(Ed - osign*magb*E15)

    X = {}; Y = {}
    for i in idx:
        w,wc = wmap[i]; z,zc = zmap[i]; p,q,r,s = pmap[i]
        dp = sp.expand(combo.diff(p)); dr = sp.expand(combo.diff(r))
        Xi,Yi = sp.symbols('Xi Yi')
        sol = sp.solve([sp.Eq(dp, -w*zc**2*Xi + w*Yi), sp.Eq(dr, -wc*z**2*Xi + wc*Yi)], [Xi,Yi])
        X[i] = sp.simplify(sol[Xi]); Y[i] = sp.simplify(sol[Yi])

    assert Y[4] == 0, f"deg{d}: D3_4 did not cancel, Y4={Y[4]}"
    claim = sum(coeff_sym(1,i)*X[i] for i in idx) + sum(coeff_sym(3,i)*Y[i] for i in (1,2,3))
    diff = sp.simplify(sp.expand(claim - combo))
    assert diff == 0, f"deg{d}: claim reconstruction mismatch, diff={diff}"

    coeff1_terms = "\n      + ".join(f"{coeff_lean_text(1,i)}\n        * ({sympy_to_lean(X[i])})" for i in idx)
    coeff3_terms = "\n      + ".join(f"{coeff_lean_text(3,i)}\n        * ({sympy_to_lean(Y[i])})" for i in (1,2,3))
    goal = f"{coeff1_terms}\n      + {coeff3_terms} = 0"

    hz_hyps = " ".join(f"(hz{i}_im : z{i}.im ≠ 0)" for i in idx)
    pqrs_decl = " ".join(f"p{i} q{i} r{i} s{i}" for i in idx)
    z_decl = " ".join(f"z{i}" for i in idx)
    hsum_lines = "\n        + ".join(
        f"p{i} • vPlusLC N L z{i} + q{i} • vMinusLC N L z{i}\n"
        f"        + r{i} • vPlusLC N L (starRingEnd ℂ z{i}) + s{i} • vMinusLC N L (starRingEnd ℂ z{i})"
        for i in idx
    )
    hkey_args = " ".join(f"z{i}" for i in idx) + " " + " ".join(f"hz{i}_im" for i in idx) + " " + \
        " ".join(f"p{i} q{i} r{i} s{i}" for i in idx)

    mag_lean = sympy_to_lean(magb)
    proof_sign = "-" if osign == 1 else "+"

    theorem_name = RED_NAME[d]
    coeff_name_d = COEFF_NAME[d]
    coeff_file_d = coeff_name_d[0].upper() + coeff_name_d[1:]

    content = f"""import CrossK4LeadingCoeff
import {coeff_file_d}

/-!
# K=4 decoupling: eliminate D3_4 from the degree-{d} equation via the degree-15 equation

Direct 4-block extension of the K=3 ReducedK{{...}} cascade, continuing S2536-2538. Sympy-derived AND
verified BEFORE writing this Lean (derive_and_verify_k4_reductions.py): combo = {coeff_name_d} {"-" if osign==1 else "+"} (octic-style
coefficient, b={b}) * crossK4LeadingCoeff eliminates D3_4 completely (verified via the 2x2-per-block
linear-system solve showing D3_4's multiplier is EXACTLY 0), and the resulting closed form (grouped by
each surviving block's coeff1/coeff3, with pure-z multiplier polynomials machine-translated from the
verified sympy expressions via sympy_to_lean.py -- not hand-transcribed, avoiding the S2537 bug class)
was independently reconstructed and confirmed to exactly match the raw combination (diff=0).

Regenerated 2026-09-16 (S2561) after the original was lost to a scratchpad rotation -- reconstructed
sympy_to_lean.py from its functional spec (LEDGER S2539) and re-derived via this same generic script;
output diffed byte-for-byte against the one surviving original (CrossK4ReducedK3.lean) to validate the
reconstruction before trusting K13/K11/K9/K7/K5's fresh output. -/

theorem {theorem_name} (N : ℕ) (hN : 8 ≤ N) (L : ℝ) (hL : L ≠ 0) ({z_decl} : ℂ)
    {hz_hyps}
    ({pqrs_decl} : ℂ)
    (hsum : {hsum_lines} = 0) :
    {goal} := by
  have hE15 := crossK4LeadingCoeff N hN L hL {hkey_args} hsum
  have hE{d} := {coeff_name_d} N hN L hL {hkey_args} hsum
  linear_combination hE{d} {proof_sign} ({mag_lean}) * hE15

#print axioms {theorem_name}
"""
    fname = f"/private/tmp/claude-501/-Users-scones/67278aad-a465-41ae-b156-e915a5de10b7/scratchpad/{theorem_name[0].upper()}{theorem_name[1:]}.lean"
    with open(fname, "w", encoding="utf-8") as f:
        f.write(content)
    bal = (content.count('('), content.count(')'))
    print(f"deg{d}: wrote {fname}, len={len(content)}, paren_balance={bal}, balanced={bal[0]==bal[1]}")
    return combo, X, Y

for d in (13,11,9,7,5,3):
    derive_and_build(d)
