import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Tactic

set_option maxHeartbeats 1000000

/-!
# §1703a — EXPONENTIAL ABSORPTION: `(a+bu)e^{-εu} ≤ a + b/ε` and `u²e^{-εu} ≤ 4/ε²`

The two absorption facts that turn the optimization's polynomial-in-`√log x` prefactors into
constants: linear via `e^v ≥ v`, quadratic via `e^v = (e^{v/2})² ≥ (1+v/2)² ≥ v²/4`.
-/

theorem linear_mul_exp_neg_le {a b ε u : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) (hε : 0 < ε)
    (hu : 0 ≤ u) : (a + b * u) * Real.exp (-(ε * u)) ≤ a + b / ε := by
  have hexp : Real.exp (-(ε * u)) ≤ 1 := by
    rw [Real.exp_le_one_iff]
    nlinarith
  have hexppos : 0 < Real.exp (-(ε * u)) := Real.exp_pos _
  rcases eq_or_lt_of_le hu with rfl | hupos
  · simp only [mul_zero, add_zero, neg_zero, Real.exp_zero, mul_one]
    have : 0 ≤ b / ε := by positivity
    linarith
  · have hkey : ε * u ≤ Real.exp (ε * u) := by
      have := Real.add_one_le_exp (ε * u)
      linarith
    have hεu : 0 < ε * u := by positivity
    have hkey2 : Real.exp (-(ε * u)) ≤ (ε * u)⁻¹ := by
      rw [Real.exp_neg, inv_eq_one_div, inv_eq_one_div]
      exact one_div_le_one_div_of_le hεu hkey
    have hbu : b * u * Real.exp (-(ε * u)) ≤ b / ε := by
      calc b * u * Real.exp (-(ε * u)) ≤ b * u * (ε * u)⁻¹ :=
            mul_le_mul_of_nonneg_left hkey2 (by positivity)
        _ = b / ε := by field_simp
    calc (a + b * u) * Real.exp (-(ε * u))
        = a * Real.exp (-(ε * u)) + b * u * Real.exp (-(ε * u)) := by ring
      _ ≤ a * 1 + b / ε := by
          apply add_le_add _ hbu
          exact mul_le_mul_of_nonneg_left hexp ha
      _ = a + b / ε := by ring

theorem sq_mul_exp_neg_le {ε u : ℝ} (hε : 0 < ε) (hu : 0 ≤ u) :
    u ^ 2 * Real.exp (-(ε * u)) ≤ 4 / ε ^ 2 := by
  have hquad : (ε * u) ^ 2 / 4 ≤ Real.exp (ε * u) := by
    have h1 : Real.exp (ε * u) = Real.exp (ε * u / 2) ^ 2 := by
      rw [sq, ← Real.exp_add]
      congr 1
      ring
    have h2 : 1 + ε * u / 2 ≤ Real.exp (ε * u / 2) := by
      have := Real.add_one_le_exp (ε * u / 2)
      linarith
    have h3 : (ε * u / 2) ^ 2 ≤ (1 + ε * u / 2) ^ 2 := by
      have hnn : 0 ≤ ε * u / 2 := by positivity
      nlinarith
    calc (ε * u) ^ 2 / 4 = (ε * u / 2) ^ 2 := by ring
      _ ≤ (1 + ε * u / 2) ^ 2 := h3
      _ ≤ Real.exp (ε * u / 2) ^ 2 := by
          have hnn : 0 ≤ 1 + ε * u / 2 := by positivity
          nlinarith [h2]
      _ = Real.exp (ε * u) := h1.symm
  rw [Real.exp_neg]
  rw [mul_inv_le_iff₀ (Real.exp_pos _), div_mul_eq_mul_div, le_div_iff₀ (by positivity)]
  calc u ^ 2 * ε ^ 2 = (ε * u) ^ 2 := by ring
    _ ≤ 4 * Real.exp (ε * u) := by linarith

#print axioms linear_mul_exp_neg_le
#print axioms sq_mul_exp_neg_le
