import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.DerivHyp
open Real
/-- Weil/archimedean building block: `1 + x²/2 ≤ cosh x` for all real x
    (cosh x − 1 = 2·sinh²(x/2) ≥ 2·(x/2)² = x²/2). -/
theorem one_add_half_sq_le_cosh (x : ℝ) : 1 + x ^ 2 / 2 ≤ Real.cosh x := by
  have hid : Real.cosh x = 2 * Real.sinh (x / 2) ^ 2 + 1 := by
    have h2 : (2 : ℝ) * (x / 2) = x := by ring
    have hc := Real.cosh_two_mul (x / 2)
    rw [h2] at hc
    have hsq := Real.cosh_sq (x / 2)
    rw [hc, hsq]; ring
  have hsinh : (x / 2) ^ 2 ≤ Real.sinh (x / 2) ^ 2 := by
    rcases le_total 0 (x / 2) with h | h
    · have hle := (Real.self_le_sinh_iff).mpr h
      nlinarith [hle, h]
    · have hle := (Real.sinh_le_self_iff).mpr h
      nlinarith [hle, h]
  rw [hid]; nlinarith [hsinh]

#print axioms one_add_half_sq_le_cosh
