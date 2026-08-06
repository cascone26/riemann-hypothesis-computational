import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Data.Real.StarOrdered
open Matrix
/-- **Machine-checked organ-2 Weil-positivity certificate**, basis φ_a(x)=e^{-x²/2}cos(ax),
a∈{15,25,35,45}.  The rigorous Guinand–Weil form Q=arch+pole−prime (operators/organ2_certify.py,
min-eig≈0.048), rationalized (rounding ≪ margin), is `PosSemidef` — verified by its exact LDL^T
sum-of-squares.  The Lean end of the Arb → exact-rational Weil certificate (§731/§761/§762). -/
theorem organ2_weil_cert_posSemidef :
    (Matrix.of ![![(14859/10000 : ℝ), (0 : ℝ), (0 : ℝ), (0 : ℝ)], ![(0 : ℝ), (7853/2500 : ℝ), (0 : ℝ), (0 : ℝ)], ![(0 : ℝ), (0 : ℝ), (481/10000 : ℝ), (0 : ℝ)], ![(0 : ℝ), (0 : ℝ), (0 : ℝ), (1917/10000 : ℝ)]]).PosSemidef := by
  rw [posSemidef_iff_dotProduct_mulVec]
  refine ⟨?_, ?_⟩
  · ext i j; fin_cases i <;> fin_cases j <;> norm_num [Matrix.conjTranspose_apply, Matrix.of_apply]
  · intro x
    simp only [dotProduct, Matrix.mulVec, Matrix.of_apply, Fin.sum_univ_four,
      Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons, Matrix.cons_val_two,
      Matrix.cons_val_three, Matrix.tail_cons, star_trivial]
    nlinarith [mul_nonneg (by norm_num : (0:ℝ) ≤ (14859/10000 : ℝ)) (sq_nonneg (x 0)),
               mul_nonneg (by norm_num : (0:ℝ) ≤ (7853/2500 : ℝ)) (sq_nonneg (x 1)),
               mul_nonneg (by norm_num : (0:ℝ) ≤ (481/10000 : ℝ)) (sq_nonneg (x 2)),
               mul_nonneg (by norm_num : (0:ℝ) ≤ (1917/10000 : ℝ)) (sq_nonneg (x 3))]

#print axioms organ2_weil_cert_posSemidef
