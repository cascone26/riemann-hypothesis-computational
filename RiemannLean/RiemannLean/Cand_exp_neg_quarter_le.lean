import Mathlib.Analysis.SpecialFunctions.Exponential
/-- Interval-lifting building block: `exp(-1/4) ≤ 4/5`, a rigorous rational bound on the Gaussian
factor `e^{-y²/4}` (here at `y²=1`) in the Weil prime term.  From `5/4 ≤ exp(1/4)` (`add_one_le_exp`)
by taking inverses.  A first step toward Lean-verified bounds on the transcendental Weil entries (§764). -/
theorem exp_neg_quarter_le : Real.exp (-(1/4)) ≤ 4/5 := by
  have h : (5:ℝ)/4 ≤ Real.exp (1/4) := by
    have := Real.add_one_le_exp (1/4 : ℝ); linarith
  rw [Real.exp_neg]
  rw [inv_le_comm₀ (Real.exp_pos _) (by norm_num)]
  linarith

#print axioms exp_neg_quarter_le
