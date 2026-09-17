import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Sqrt

/-!
# §1445 — threshold `C·log x ≤ x` for `x ≥ (2C)²` (`C ≥ 0`)

`log x = 2·log√x ≤ 2(√x - 1) < 2√x`, and `x ≥ (2C)²` gives `√x ≥ 2C`, so `2C·√x ≤ √x·√x = x`.  Hence
`C·log x ≤ 2C√x ≤ x`.  This is the scalar threshold that closes the constant-choice fixed-point: with `x = |γ|`,
`C = 512`, once `|γ| ≥ (2·512)²` the self-referential BC inequality `B ≥ …log|γ|…` is satisfiable.
EXIT0 / std-axioms.
-/

theorem mul_log_le_self {C x : ℝ} (hC : 0 ≤ C) (hx : (2 * C) ^ 2 ≤ x) : C * Real.log x ≤ x := by
  have hx0 : (0 : ℝ) ≤ x := le_trans (by positivity) hx
  rcases hx0.eq_or_lt with h | hxpos
  · rw [← h]; simp
  · have hlogsq : Real.log x = 2 * Real.log (Real.sqrt x) := by
      rw [Real.log_sqrt hx0]; ring
    have hsqrt : Real.log x ≤ 2 * (Real.sqrt x - 1) := by
      have hle := Real.log_le_sub_one_of_pos (Real.sqrt_pos.mpr hxpos)
      rw [hlogsq]; linarith
    have hge : 2 * C ≤ Real.sqrt x := by
      rw [show (2 * C : ℝ) = Real.sqrt ((2 * C) ^ 2) from (Real.sqrt_sq (by positivity)).symm]
      exact Real.sqrt_le_sqrt hx
    have h2Csqrt : 2 * C * Real.sqrt x ≤ x := by
      calc 2 * C * Real.sqrt x ≤ Real.sqrt x * Real.sqrt x := by
            nlinarith [Real.sqrt_nonneg x, hge]
        _ = x := Real.mul_self_sqrt hx0
    nlinarith [mul_le_mul_of_nonneg_left hsqrt hC, h2Csqrt, hC]

#print axioms mul_log_le_self
