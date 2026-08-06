import Mathlib.Analysis.SpecialFunctions.Exponential
/-- Interval-lifting brick: Gaussian tail bound `e^{−y²/4} ≤ e^{−y/4}` for `y ≥ 1`
(since `y² ≥ y` there).  Gives a summable bound for the Weil prime-sum tail (LEDGER §732/§765). -/
theorem gaussian_le_exp_neg_quarter {y : ℝ} (hy : 1 ≤ y) :
    Real.exp (-(y ^ 2) / 4) ≤ Real.exp (-y / 4) := by
  apply Real.exp_le_exp.mpr
  have hyy : y ≤ y ^ 2 := by nlinarith [hy]
  linarith

#print axioms gaussian_le_exp_neg_quarter
