import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
open Real
/-- A fresh tractable Handle: for 0 ≤ x, cosh x ≥ 1 + x^2/2 is false in general? No—
    prove the clean true fact: cos x ≤ 1 and hence 2 - 2*cos x ≥ 0 (a Weil-form-flavored kernel). -/
theorem two_sub_two_cos_nonneg (x : ℝ) : 0 ≤ 2 - 2 * cos x := by
  have := Real.cos_le_one x
  linarith
