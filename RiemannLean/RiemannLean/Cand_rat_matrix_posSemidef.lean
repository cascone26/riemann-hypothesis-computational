import Mathlib.LinearAlgebra.Matrix.PosDef
open Matrix
/-- Concept proof for the Lean-verified certificate bridge (LEDGER §755/§758): a concrete PSD
matrix with rational entries is machine-checked `PosSemidef` over ℝ (quadratic form = sum of
squares) — the Lean end of a rigorous finite Weil-positivity certificate. -/
theorem rat_matrix_posSemidef :
    (Matrix.of ![![(2 : ℝ), 1], ![1, 2]]).PosSemidef := by
  rw [posSemidef_iff_dotProduct_mulVec]
  refine ⟨?_, ?_⟩
  · ext i j; fin_cases i <;> fin_cases j <;> rfl
  · intro x
    simp only [dotProduct, Matrix.mulVec, Matrix.of_apply, Fin.sum_univ_two,
      Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons, star_trivial]
    nlinarith [sq_nonneg (x 0 + x 1), sq_nonneg (x 0), sq_nonneg (x 1)]

#print axioms rat_matrix_posSemidef
