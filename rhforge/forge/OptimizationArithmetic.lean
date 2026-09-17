import ExpAbsorption
import Mathlib.Analysis.SpecialFunctions.Pow.Real

set_option maxHeartbeats 1000000

/-!
# §1703c — OPTIMIZATION ARITHMETIC: the rpow/log/absorption toolkit for `σ₀ = 1+1/log x`, `T = e^{√log x}`

Every arithmetic identity the §1703 instantiation needs, isolated: `x^{σ₀+1} = e·x²`, the
`σ₁`-power collapse `x^{σ₁+1} ≤ e^{2c₁}·x²·e^{-c₁u}`, the `log(e^u+2) ∈ [u, u+2]` bracket, and
the three absorption combos that turn `(3+u)`, `(u+2)²`, `(u²+A)` prefactors into constants.
-/

theorem rpow_sigma0_add_one_eq {x : ℝ} (hx1 : 1 < x) :
    x ^ (1 + 1 / Real.log x + 1) = Real.exp 1 * x ^ (2 : ℝ) := by
  have hx0 : (0 : ℝ) < x := by linarith
  have hlog : 0 < Real.log x := Real.log_pos hx1
  have h1 : 1 + 1 / Real.log x + 1 = (2 : ℝ) + 1 / Real.log x := by ring
  rw [h1, Real.rpow_add hx0, Real.rpow_def_of_pos hx0 (1 / Real.log x),
    mul_one_div, div_self hlog.ne']
  ring

theorem rpow_sigma1_add_one_le {x T c₁ : ℝ} (hc₁ : 0 < c₁) (hx1 : 1 < x)
    (hlog1 : 1 ≤ Real.log x)
    (hLpos : 0 < Real.log (T + 2))
    (hLhi : Real.log (T + 2) ≤ Real.sqrt (Real.log x) + 2) :
    x ^ (1 - c₁ / Real.log (T + 2) + 1)
      ≤ Real.exp (2 * c₁) * x ^ (2 : ℝ)
          * Real.exp (-(c₁ * Real.sqrt (Real.log x))) := by
  have hx0 : (0 : ℝ) < x := by linarith
  have hu0 : 0 ≤ Real.sqrt (Real.log x) := Real.sqrt_nonneg _
  have hu2 : Real.sqrt (Real.log x) ^ 2 = Real.log x := Real.sq_sqrt (by linarith)
  set u := Real.sqrt (Real.log x) with hudef
  set L := Real.log (T + 2) with hLdef
  have h1 : 1 - c₁ / L + 1 = (2 : ℝ) + (-(c₁ / L)) := by ring
  rw [h1, Real.rpow_add hx0, Real.rpow_def_of_pos hx0 (-(c₁ / L))]
  have hkey : Real.log x * (-(c₁ / L)) ≤ 2 * c₁ + (-(c₁ * u)) := by
    rw [← hu2]
    have h2 : (c₁ * u - 2 * c₁) * L ≤ c₁ * u ^ 2 := by
      rcases le_total u 2 with h | h
      · have hnp : c₁ * u - 2 * c₁ ≤ 0 := by nlinarith
        nlinarith [hnp, hLpos.le, mul_nonneg hc₁.le (sq_nonneg u)]
      · have hprod : (0 : ℝ) ≤ (u - 2) * (u + 2 - L) :=
          mul_nonneg (by linarith) (by linarith)
        nlinarith [mul_nonneg hc₁.le hprod]
    have h3 : u ^ 2 * (-(c₁ / L)) = -(c₁ * u ^ 2 / L) := by ring
    rw [h3, neg_le, le_div_iff₀ hLpos]
    nlinarith [h2]
  calc x ^ (2 : ℝ) * Real.exp (Real.log x * (-(c₁ / L)))
      ≤ x ^ (2 : ℝ) * Real.exp (2 * c₁ + (-(c₁ * u))) :=
        mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr hkey) (Real.rpow_nonneg hx0.le _)
    _ = Real.exp (2 * c₁) * x ^ (2 : ℝ) * Real.exp (-(c₁ * u)) := by
        rw [Real.exp_add]; ring

theorem log_exp_add_two_bounds {u : ℝ} (hu : 0 ≤ u) :
    u ≤ Real.log (Real.exp u + 2) ∧ Real.log (Real.exp u + 2) ≤ u + 2 := by
  have hepos : (0 : ℝ) < Real.exp u := Real.exp_pos u
  constructor
  · calc u = Real.log (Real.exp u) := (Real.log_exp u).symm
      _ ≤ Real.log (Real.exp u + 2) := Real.log_le_log hepos (by linarith)
  · have h3 : Real.exp u + 2 ≤ 3 * Real.exp u := by
      have h1 : (1 : ℝ) ≤ Real.exp u := by
        have := Real.add_one_le_exp u
        linarith
      linarith
    have hlog3 : Real.log 3 ≤ 2 := by
      have h3e : (3 : ℝ) ≤ Real.exp 2 := by
        have := Real.add_one_le_exp (2 : ℝ)
        linarith
      have := Real.log_le_log (by norm_num : (0 : ℝ) < 3) h3e
      rwa [Real.log_exp] at this
    calc Real.log (Real.exp u + 2)
        ≤ Real.log (3 * Real.exp u) := Real.log_le_log (by linarith) h3
      _ = Real.log 3 + u := by
          rw [Real.log_mul (by norm_num) hepos.ne', Real.log_exp]
      _ ≤ u + 2 := by linarith

theorem linear_absorb_half {c₁ u : ℝ} (hc₁ : 0 < c₁) (hu : 0 ≤ u) :
    (3 + u) * Real.exp (-(c₁ * u)) ≤ (3 + 2 / c₁) * Real.exp (-(c₁ / 2 * u)) := by
  have hsplit : Real.exp (-(c₁ * u))
      = Real.exp (-(c₁ / 2 * u)) * Real.exp (-(c₁ / 2 * u)) := by
    rw [← Real.exp_add]; congr 1; ring
  rw [hsplit, ← mul_assoc]
  apply mul_le_mul_of_nonneg_right _ (Real.exp_pos _).le
  have h := linear_mul_exp_neg_le (a := 3) (b := 1) (ε := c₁ / 2) (u := u)
    (by norm_num) (by norm_num) (by positivity) hu
  calc (3 + u) * Real.exp (-(c₁ / 2 * u))
      = (3 + 1 * u) * Real.exp (-(c₁ / 2 * u)) := by ring
    _ ≤ 3 + 1 / (c₁ / 2) := h
    _ = 3 + 2 / c₁ := by rw [one_div_div]

theorem sq_shift_absorb {u : ℝ} (hu : 0 ≤ u) :
    (u + 2) ^ 2 * Real.exp (-(2 * u)) ≤ 12 * Real.exp (-u) := by
  have hsplit : Real.exp (-(2 * u)) = Real.exp (-u) * Real.exp (-u) := by
    rw [← Real.exp_add]; congr 1; ring
  rw [hsplit, ← mul_assoc]
  apply mul_le_mul_of_nonneg_right _ (Real.exp_pos _).le
  have h1 := sq_mul_exp_neg_le (ε := 1) (u := u) one_pos hu
  have h2 := linear_mul_exp_neg_le (a := 4) (b := 4) (ε := 1) (u := u)
    (by norm_num) (by norm_num) one_pos hu
  simp only [one_mul] at h1 h2
  norm_num at h1 h2
  nlinarith [h1, h2, Real.exp_pos (-u)]

theorem sq_plus_const_absorb {A u : ℝ} (hA : 0 ≤ A) (hu : 0 ≤ u) :
    (u ^ 2 + A) * Real.exp (-u) ≤ (16 + A) * Real.exp (-(u / 2)) := by
  have hsplit : Real.exp (-u) = Real.exp (-(u / 2)) * Real.exp (-(u / 2)) := by
    rw [← Real.exp_add]; congr 1; ring
  rw [hsplit, ← mul_assoc]
  apply mul_le_mul_of_nonneg_right _ (Real.exp_pos _).le
  have h1 := sq_mul_exp_neg_le (ε := 1 / 2) (u := u) (by norm_num) hu
  have he : -(1 / 2 * u) = -(u / 2) := by ring
  rw [he] at h1
  norm_num at h1
  have hexp1 : Real.exp (-(u / 2)) ≤ 1 := by
    rw [Real.exp_le_one_iff]
    linarith
  have hAE : A * Real.exp (-(u / 2)) ≤ A := mul_le_of_le_one_right hA hexp1
  nlinarith [h1, hAE, Real.exp_pos (-(u / 2))]

#print axioms rpow_sigma0_add_one_eq
#print axioms rpow_sigma1_add_one_le
#print axioms log_exp_add_two_bounds
#print axioms linear_absorb_half
#print axioms sq_shift_absorb
#print axioms sq_plus_const_absorb
