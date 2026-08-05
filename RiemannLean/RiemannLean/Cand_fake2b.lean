import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
theorem fake2 (x : ℝ) : Real.cos x ≤ 0 := by
  have := Real.cos_le_one x
  linarith

#print axioms fake2
