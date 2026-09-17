import Mathlib.Analysis.SpecialFunctions.Exp

/-!
# §1293 — `e^y - 1 ≤ y·e^y` for `y ≥ 0` (the divided-difference collapse)

The convexity inequality behind turning §1291's first main term `(N^{1-σ}-1)/(1-σ)` into `≤ N^{1-σ}·log N`.
Writing `y = (1-σ)·log N ≥ 0`, the divided difference `(N^{1-σ}-1)/(1-σ) = ((e^y-1)/y)·log N`, and `(e^y-1)/y ≤ e^y`
is exactly `e^y - 1 ≤ y·e^y`.  From `Real.add_one_le_exp (-y)` (`1-y ≤ e^{-y}`) times `e^y`.  EXIT0 / std-axioms.
-/

theorem exp_sub_one_le_mul_exp {y : ℝ} (_hy : 0 ≤ y) : Real.exp y - 1 ≤ y * Real.exp y := by
  have h : -y + 1 ≤ Real.exp (-y) := Real.add_one_le_exp (-y)
  have hey : (0 : ℝ) < Real.exp y := Real.exp_pos y
  have key : (-y + 1) * Real.exp y ≤ Real.exp (-y) * Real.exp y :=
    mul_le_mul_of_nonneg_right h (le_of_lt hey)
  rw [← Real.exp_add, neg_add_cancel, Real.exp_zero] at key
  nlinarith [key]

#print axioms exp_sub_one_le_mul_exp
