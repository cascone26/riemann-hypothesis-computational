import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Data.Real.StarOrdered
open Matrix
/-- **Concrete Weil-style positivity certificate (3×3).** A rational SPD matrix is `PosSemidef`,
proven via its LDL^T sum-of-squares: `xᵀMx = (3/2)(x₀+x₁/6)² + (59/24)(x₁+8x₂/59)² + (346/177)x₂²`.
This is the concrete Lean end of an organ-2 Arb → exact-rational-LDL^T Weil certificate. -/
theorem weil_style_cert_3x3 :
    (Matrix.of ![![(3/2 : ℝ), 1/4, 0], ![1/4, 5/2, 1/3], ![0, 1/3, 2]]).PosSemidef := by
  rw [posSemidef_iff_dotProduct_mulVec]
  refine ⟨?_, ?_⟩
  · ext i j; fin_cases i <;> fin_cases j <;> norm_num [Matrix.conjTranspose_apply, Matrix.of_apply]
  · intro x
    simp only [dotProduct, Matrix.mulVec, Matrix.of_apply, Fin.sum_univ_three,
      Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons, Matrix.cons_val_two,
      Matrix.tail_cons, star_trivial]
    nlinarith [mul_nonneg (by norm_num : (0:ℝ) ≤ 3/2) (sq_nonneg (x 0 + x 1 / 6)),
               mul_nonneg (by norm_num : (0:ℝ) ≤ 59/24) (sq_nonneg (x 1 + 8 * x 2 / 59)),
               mul_nonneg (by norm_num : (0:ℝ) ≤ 346/177) (sq_nonneg (x 2)),
               sq_nonneg (x 0), sq_nonneg (x 1), sq_nonneg (x 2)]

#print axioms weil_style_cert_3x3
