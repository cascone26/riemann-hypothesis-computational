#!/usr/bin/env python3
"""Lean replay of the ZERO-FREE open-rung certificate (§788 → §792; lane (e) per §788-addendum:
"the Lean LDLᵀ rational replay should replay THESE unconditional certificates").

Recomputes the 4x4 zero-free Weil Gram at 4W=1.32, span {14,22,30,38}, RC=400 (pole + arch − prime,
Arb balls — the §788 route, arb_zerofree_cert3 machinery), rationalizes the certified midpoints,
bounds ‖Q−M‖₂ ≤ ρ (radii + rationalization, Frobenius), takes exact LDLᵀ of M over ℚ, and emits a
Lean theorem proving M.PosSemidef by sum-of-squares — so the kernel-checked object carries the
HONEST hypothesis set (none: no zero knowledge in this route).  Pattern: gen_lean_weil_cert_rigorous
(§767), ported from the Gaussian organ-2 basis to the Fejér zero-free route.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import sympy as sp
from sympy import Rational
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc

fourW = 1.32; CENTERS = [14, 22, 30, 38]; RC = 400
ctx.prec = 220
ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
cs = [arb(c) for c in CENTERS]; nb = len(cs)

print(f"computing zero-free Arb Gram (4W={fourW}, span {CENTERS}, RC={RC})...", flush=True)
Qb = [[None]*nb for _ in range(nb)]
maxrad = arb(0)
for j in range(nb):
    for k in range(j, nb):
        v = ac.pole_entry(cs[j], cs[k]) + zc.arch_entry_v3(cs[j], cs[k], arb(RC)) - ac.prime_entry(cs[j], cs[k])
        Qb[j][k] = Qb[k][j] = v
        if v.rad() > maxrad: maxrad = v.rad()
    print(f"  row {j} done (diag rad {float(Qb[j][j].rad()):.2e})", flush=True)

M = sp.zeros(nb, nb); radsq = arb(0)
for j in range(nb):
    for k in range(nb):
        mid = Qb[j][k].mid(); rad = Qb[j][k].rad()
        M[j, k] = Rational(str(float(mid))).limit_denominator(10**6)
        re = abs(float(mid) - float(M[j, k]))
        radsq += arb(re + float(rad))**2
rho = radsq.sqrt(); rho_up = float(rho) + 1e-15
print(f"  ARB max entry radius: {float(maxrad):.2e}")
print(f"  certified ||Q-M||_2 <= rho <= {rho_up:.3e}")
L, D = M.LDLdecomposition(); d = [D[i, i] for i in range(nb)]
print(f"  M pivots (exact): min = {min(float(x) for x in d):.5f}  all>0: {all(x > 0 for x in d)}")
# rigorous min-eig(M) > rho: exact LDL^T of (M - rho_Q*I) over Q — all pivots > 0 proves it,
# no float eigensolver anywhere (stronger than the §767 generator's float min-eig print).
rho_Q = Rational(rho_up).limit_denominator(10**12)
if float(rho_Q) < rho_up: rho_Q = rho_Q + Rational(1, 10**9)
Ms = M - rho_Q * sp.eye(nb)
Ls, Ds = Ms.LDLdecomposition(); ds = [Ds[i, i] for i in range(nb)]
assert all(x > 0 for x in d) and all(x > 0 for x in ds)
mineigM = min(float(x) for x in d)   # display only; the rigorous statement is the (M - rho I) LDL^T
print(f"  (M - rho*I) exact pivots all > 0: {all(x > 0 for x in ds)} => min-eig(M) > rho => Q PD RIGOROUSLY")

LT = L.T
def lr(x):
    x = sp.Rational(x); return f"({x.p}/{x.q} : ℝ)" if x.q != 1 else f"({x.p} : ℝ)"
def ell(i):
    ts = [("x " + str(j) if LT[i, j] == 1 else f"{lr(LT[i, j])} * x {j}") for j in range(nb) if LT[i, j] != 0]
    return " + ".join(ts) if ts else "0"
def row(i): return "![" + ", ".join(lr(M[i, j]) for j in range(nb)) + "]"
hints = ",\n               ".join(
    f"mul_nonneg (by norm_num : (0:ℝ) ≤ {lr(d[i])}) (sq_nonneg ({ell(i)}))" for i in range(nb))
lean = f'''import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Data.Real.StarOrdered
open Matrix
/-- **Machine-checked ZERO-FREE open-rung Weil certificate matrix** (Fejér span c∈{{14,22,30,38}},
support 4W=1.32 > log 3 — the first open rung, primes {{2,3}} in support).
`M` rationalizes the ARB-CERTIFIED midpoints of the zero-free route entries pole+arch−prime
(operators/arb_zerofree_cert3.py, §788: NO zero knowledge anywhere — no verified-zero list, no
on-line-ness hypothesis, no strip lemma). This theorem proves `M.PosSemidef` via its exact LDLᵀ
sum-of-squares. With the certified bound ‖Q−M‖₂ ≤ ρ ≈ {rho_up:.2e} and min-eig(M) > ρ (proven by
exact LDLᵀ of M − ρI over ℚ, all pivots positive; min LDLᵀ pivot of M = {mineigM:.4f}), the TRUE
zero-free Weil form is positive-definite on this span UNCONDITIONALLY (§788/§792). -/
theorem zerofree_rung_certM_posSemidef :
    (Matrix.of ![{row(0)}, {row(1)}, {row(2)}, {row(3)}]).PosSemidef := by
  rw [posSemidef_iff_dotProduct_mulVec]
  refine ⟨?_, ?_⟩
  · ext i j; fin_cases i <;> fin_cases j <;> norm_num [Matrix.conjTranspose_apply, Matrix.of_apply]
  · intro x
    simp only [dotProduct, Matrix.mulVec, Matrix.of_apply, Fin.sum_univ_four,
      Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons, Matrix.cons_val_two,
      Matrix.cons_val_three, Matrix.tail_cons, star_trivial]
    nlinarith [{hints}]
'''
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "RiemannLean", "RiemannLean",
                   "Cand_zerofree_rung_certM.lean")
open(out, "w").write(lean)
print(f"emitted {out}")
