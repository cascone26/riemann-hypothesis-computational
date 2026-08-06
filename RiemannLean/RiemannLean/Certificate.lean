import Mathlib.LinearAlgebra.Matrix.PosDef

/-!
# The Lean–Arb certificate bridge (concept)

LEDGER §755 identified the boundary: the Weil form's archimedean term is a *transcendental*
integral, so pure Lean cannot prove positivity symbolically — it needs *interval-arithmetic*
certification (Arb, organ 2). But the two lanes DO compose: Arb rigorously encloses the
transcendental Weil-form entries in a rational lower-bound matrix `M` and certifies `M ⪰ 0` by
exact-rational LDL^T; **Lean then machine-checks that `M` is `PosSemidef`** (its quadratic form is
a sum of squares). The one step Lean cannot yet check is the enclosure `M ⪯ (true Weil matrix)`,
which is the Arb (transcendental) part — a known gap requiring a Lean-verified interval layer.

This module demonstrates the Lean end: a concrete rational-entry PSD matrix is machine-checked
`PosSemidef` over ℝ via `posSemidef_iff_dotProduct_mulVec` + a sum-of-squares witness. This is the
prototype for lifting an organ-2 rational Weil certificate (`operators/organ2_certify.py`,
`weil_cert_15_25_35_45.json`) into a machine-checked Lean proof.
-/

open Matrix

/-- **Certificate-bridge concept proof.**  A concrete rational-entry matrix is `PosSemidef` over ℝ,
machine-verified: its quadratic form `2x₀² + 2x₀x₁ + 2x₁² = (x₀+x₁)² + x₀² + x₁²` is a sum of
squares.  This is the Lean end of a rigorous finite Weil-positivity certificate. -/
theorem rat_matrix_posSemidef :
    (Matrix.of ![![(2 : ℝ), 1], ![1, 2]]).PosSemidef := by
  rw [posSemidef_iff_dotProduct_mulVec]
  refine ⟨?_, ?_⟩
  · ext i j; fin_cases i <;> fin_cases j <;> rfl
  · intro x
    simp only [dotProduct, Matrix.mulVec, Matrix.of_apply, Fin.sum_univ_two,
      Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons, star_trivial]
    nlinarith [sq_nonneg (x 0 + x 1), sq_nonneg (x 0), sq_nonneg (x 1)]
