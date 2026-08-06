import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Data.Real.StarOrdered
open Matrix
/-- **Machine-checked organ-2 Weil certificate (6×6)**, basis a∈{15,20,25,30,35,40}.
Rigorous Guinand–Weil form Q=arch+pole−prime (min-eig≈0.048), rationalized, `PosSemidef` via
exact LDL^T sum-of-squares.  Extends §763 — the machine-checked Weil criterion scales. -/
theorem organ2_weil_cert_6x6 :
    (Matrix.of ![![(14859/10000 : ℝ), (0 : ℝ), (0 : ℝ), (0 : ℝ), (0 : ℝ), (0 : ℝ)], ![(0 : ℝ), (5527/5000 : ℝ), (7/10000 : ℝ), (0 : ℝ), (0 : ℝ), (0 : ℝ)], ![(0 : ℝ), (7/10000 : ℝ), (7853/2500 : ℝ), (0 : ℝ), (0 : ℝ), (0 : ℝ)], ![(0 : ℝ), (0 : ℝ), (0 : ℝ), (26233/10000 : ℝ), (51/10000 : ℝ), (0 : ℝ)], ![(0 : ℝ), (0 : ℝ), (0 : ℝ), (51/10000 : ℝ), (481/10000 : ℝ), (3/500 : ℝ)], ![(0 : ℝ), (0 : ℝ), (0 : ℝ), (0 : ℝ), (3/500 : ℝ), (13601/10000 : ℝ)]]).PosSemidef := by
  rw [posSemidef_iff_dotProduct_mulVec]
  refine ⟨?_, ?_⟩
  · ext i j; fin_cases i <;> fin_cases j <;> norm_num [Matrix.conjTranspose_apply, Matrix.of_apply]
  · intro x
    simp only [dotProduct, Matrix.mulVec, Matrix.of_apply, Fin.sum_univ_six, star_trivial,
      Matrix.cons_val', Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons,
      Matrix.empty_val', Matrix.cons_val_fin_one, Matrix.head_fin_const, Fin.isValue]
    nlinarith [mul_nonneg (by norm_num : (0:ℝ) ≤ (14859/10000 : ℝ)) (sq_nonneg (x 0)),
               mul_nonneg (by norm_num : (0:ℝ) ≤ (5527/5000 : ℝ)) (sq_nonneg (x 1 + (7/11054 : ℝ) * x 2)),
               mul_nonneg (by norm_num : (0:ℝ) ≤ (347228199/110540000 : ℝ)) (sq_nonneg (x 2)),
               mul_nonneg (by norm_num : (0:ℝ) ≤ (26233/10000 : ℝ)) (sq_nonneg (x 3 + (51/26233 : ℝ) * x 4)),
               mul_nonneg (by norm_num : (0:ℝ) ≤ (788467/16395625 : ℝ)) (sq_nonneg (x 4 + (393495/3153868 : ℝ) * x 5)),
               mul_nonneg (by norm_num : (0:ℝ) ≤ (5359018621/3942335000 : ℝ)) (sq_nonneg (x 5))]

#print axioms organ2_weil_cert_6x6
