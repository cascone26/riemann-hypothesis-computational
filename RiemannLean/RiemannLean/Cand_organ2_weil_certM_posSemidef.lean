import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Data.Real.StarOrdered
open Matrix
/-- **Machine-checked organ-2 Weil certificate matrix M** (basis a∈{15,25,35,45}).
`M` is the rational matrix obtained by rationalizing the ARB-CERTIFIED midpoints of the rigorous
Guinand–Weil form entries (organ2_certify.py: each entry an Arb ball, radius ≤ 5e-21). This theorem
proves `M.PosSemidef` via its exact LDL^T sum-of-squares.  With `weil_posSemidef_of_dominates` and
the certified bound ‖Q−M‖₂ ≤ ρ ≈ 6.20e-11 ≪ min-eig(M) ≈ 0.0481, the TRUE Weil form
Q=arch+pole−prime is positive-definite (the remaining `(Q−M)⪰0` step is the Arb enclosure, §764/§766). -/
theorem organ2_weil_certM_posSemidef :
    (Matrix.of ![![(1221196/821849 : ℝ), (0 : ℝ), (0 : ℝ), (0 : ℝ)], ![(0 : ℝ), (2161912/688239 : ℝ), (0 : ℝ), (0 : ℝ)], ![(0 : ℝ), (0 : ℝ), (14289/297049 : ℝ), (0 : ℝ)], ![(0 : ℝ), (0 : ℝ), (0 : ℝ), (169601/884868 : ℝ)]]).PosSemidef := by
  rw [posSemidef_iff_dotProduct_mulVec]
  refine ⟨?_, ?_⟩
  · ext i j; fin_cases i <;> fin_cases j <;> norm_num [Matrix.conjTranspose_apply, Matrix.of_apply]
  · intro x
    simp only [dotProduct, Matrix.mulVec, Matrix.of_apply, Fin.sum_univ_four,
      Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons, Matrix.cons_val_two,
      Matrix.cons_val_three, Matrix.tail_cons, star_trivial]
    nlinarith [mul_nonneg (by norm_num : (0:ℝ) ≤ (1221196/821849 : ℝ)) (sq_nonneg (x 0)),
               mul_nonneg (by norm_num : (0:ℝ) ≤ (2161912/688239 : ℝ)) (sq_nonneg (x 1)),
               mul_nonneg (by norm_num : (0:ℝ) ≤ (14289/297049 : ℝ)) (sq_nonneg (x 2)),
               mul_nonneg (by norm_num : (0:ℝ) ≤ (169601/884868 : ℝ)) (sq_nonneg (x 3))]

#print axioms organ2_weil_certM_posSemidef
