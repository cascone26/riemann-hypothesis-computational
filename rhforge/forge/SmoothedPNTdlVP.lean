import SmoothedPerronErrorSharp
import LSeriesVonMangoldtPoleBound
import OptimizationArithmetic

set_option maxHeartbeats 4000000

open Complex MeasureTheory Chebyshev

open scoped LSeries.notation ArithmeticFunction

/-!
# §1703 ★★★★★ — DE LA VALLÉE POUSSIN'S SMOOTHED PNT: `Ψ(x) = x²/2 + O(x²·e^{-c√log x})`

The optimization pass over the free-T crown (§1702d): instantiate at `σ₀ = 1 + 1/log x`,
`T = e^{√log x}`, and absorb every polynomial-in-`√log x` prefactor via §1703a/§1703c.  With
`u = √log x`: `x^{σ₀+1} = e·x²` exactly, `log(T+2) ∈ [u, u+2]`, the `σ₁`-power collapses to
`e^{2c₁}·x²·e^{-c₁u}`, `‖L(Λ)(σ₀)‖ ≤ log x + A = u² + A` (§1703b), and the three error terms are
absorbed into `C·x²·e^{-cu}` with `c = min(c₁/2, 1/2)`.  This is the classical smoothed prime
number theorem with de la Vallée Poussin's error term, fully machine-checked from Mathlib axioms.
-/

theorem smoothed_pnt_dlvp :
    ∃ C c x₀ : ℝ, 0 < C ∧ 0 < c ∧ 1 ≤ x₀ ∧ ∀ x : ℝ, x₀ ≤ x →
      ‖(2 * Real.pi : ℝ) • psi1C x - 2 * (Real.pi : ℂ) * ((x : ℂ) ^ (2 : ℂ) / 2)‖
        ≤ C * x ^ (2 : ℝ) * Real.exp (-(c * Real.sqrt (Real.log x))) := by
  obtain ⟨c₁, K, T₂, hc₁, hK, hT₂, hmain⟩ := smoothed_perron_error_free_T_sharp
  obtain ⟨A, hA0, hA⟩ := norm_LSeries_vonMangoldt_pole_bound
  refine ⟨K * (Real.exp (2 * c₁) * (3 + 2 / c₁) + 12 * Real.exp 1 + (16 + A) * Real.exp 1),
    min (c₁ / 2) (1 / 2), Real.exp (max ((Real.log T₂) ^ 2) 1),
    by positivity, lt_min (by positivity) (by norm_num), ?_, ?_⟩
  · have h1 := Real.add_one_le_exp (max ((Real.log T₂) ^ 2) 1)
    have h2 : (0 : ℝ) ≤ max ((Real.log T₂) ^ 2) 1 :=
      le_trans zero_le_one (le_max_right _ _)
    linarith
  intro x hx
  have hx0 : (0 : ℝ) < x := lt_of_lt_of_le (Real.exp_pos _) hx
  have hlogx : max ((Real.log T₂) ^ 2) 1 ≤ Real.log x :=
    (Real.le_log_iff_exp_le hx0).mpr hx
  have hlog1 : 1 ≤ Real.log x := le_trans (le_max_right _ _) hlogx
  have hlogpos : (0 : ℝ) < Real.log x := by linarith
  have hx1 : 1 < x := by
    rcases le_total x 1 with h | h
    · have := Real.log_nonpos hx0.le h
      linarith
    · rcases eq_or_lt_of_le h with rfl | h'
      · simp at hlogpos
      · exact h'
  have hu0 : 0 ≤ Real.sqrt (Real.log x) := Real.sqrt_nonneg _
  have hu2 : Real.sqrt (Real.log x) ^ 2 = Real.log x := Real.sq_sqrt hlogpos.le
  have hu1 : 1 ≤ Real.sqrt (Real.log x) := by
    have h := Real.sqrt_le_sqrt hlog1
    rwa [Real.sqrt_one] at h
  set u := Real.sqrt (Real.log x) with hudef
  -- T = e^u qualifies
  have hT₂pos : (0 : ℝ) < T₂ := by linarith
  have hTT₂ : T₂ ≤ Real.exp u := by
    have h1 : (Real.log T₂) ^ 2 ≤ Real.log x := le_trans (le_max_left _ _) hlogx
    have h2 : Real.log T₂ ≤ u := by
      have h3 := Real.sqrt_le_sqrt h1
      rw [Real.sqrt_sq_eq_abs] at h3
      exact le_trans (le_abs_self _) h3
    calc T₂ = Real.exp (Real.log T₂) := (Real.exp_log hT₂pos).symm
      _ ≤ Real.exp u := Real.exp_le_exp.mpr h2
  -- σ₀ = 1 + 1/log x qualifies
  have hσ₀1 : 1 < 1 + 1 / Real.log x := by
    have := div_pos one_pos hlogpos
    linarith
  have hσ₀2 : 1 + 1 / Real.log x ≤ 2 := by
    have := (div_le_one hlogpos).mpr hlog1
    linarith
  have hmain' := hmain (Real.exp u) hTT₂ x hx1.le (1 + 1 / Real.log x) hσ₀1 hσ₀2
  -- log(T+2) bracket
  obtain ⟨hLlo, hLhi⟩ := log_exp_add_two_bounds hu0
  have hLpos : 0 < Real.log (Real.exp u + 2) := lt_of_lt_of_le (by linarith) hLlo
  -- x^{σ₀+1} = e·x²
  have hxσ : x ^ (1 + 1 / Real.log x + 1) = Real.exp 1 * x ^ (2 : ℝ) :=
    rpow_sigma0_add_one_eq hx1
  have hxp2 : (0 : ℝ) < x ^ (2 : ℝ) := Real.rpow_pos_of_pos hx0 _
  set c := min (c₁ / 2) (1 / 2) with hcdef
  have hc_le1 : c ≤ c₁ / 2 := min_le_left _ _
  have hc_le2 : c ≤ 1 / 2 := min_le_right _ _
  set X := x ^ (2 : ℝ) * Real.exp (-(c * u)) with hXdef
  have hX0 : (0 : ℝ) ≤ X := by positivity
  -- TERM 1: (1+log(T+2))·x^{σ₁+1} ≤ e^{2c₁}(3+2/c₁)·X
  have ht1 : (1 + Real.log (Real.exp u + 2))
      * x ^ (1 - c₁ / Real.log (Real.exp u + 2) + 1)
      ≤ Real.exp (2 * c₁) * (3 + 2 / c₁) * X := by
    have hB := rpow_sigma1_add_one_le hc₁ hx1 hlog1 hLpos hLhi
    rw [← hudef] at hB
    have hfac : 1 + Real.log (Real.exp u + 2) ≤ 3 + u := by linarith
    have hxp1 : (0 : ℝ) ≤ x ^ (1 - c₁ / Real.log (Real.exp u + 2) + 1) :=
      (Real.rpow_pos_of_pos hx0 _).le
    calc (1 + Real.log (Real.exp u + 2))
          * x ^ (1 - c₁ / Real.log (Real.exp u + 2) + 1)
        ≤ (3 + u) * (Real.exp (2 * c₁) * x ^ (2 : ℝ) * Real.exp (-(c₁ * u))) :=
          mul_le_mul hfac hB hxp1 (by linarith)
      _ = Real.exp (2 * c₁) * x ^ (2 : ℝ) * ((3 + u) * Real.exp (-(c₁ * u))) := by ring
      _ ≤ Real.exp (2 * c₁) * x ^ (2 : ℝ)
            * ((3 + 2 / c₁) * Real.exp (-(c₁ / 2 * u))) :=
          mul_le_mul_of_nonneg_left (linear_absorb_half hc₁ hu0) (by positivity)
      _ ≤ Real.exp (2 * c₁) * x ^ (2 : ℝ)
            * ((3 + 2 / c₁) * Real.exp (-(c * u))) := by
          apply mul_le_mul_of_nonneg_left _ (by positivity)
          apply mul_le_mul_of_nonneg_left _ (by positivity)
          apply Real.exp_le_exp.mpr
          have := mul_le_mul_of_nonneg_right hc_le1 hu0
          linarith
      _ = Real.exp (2 * c₁) * (3 + 2 / c₁) * X := by rw [hXdef]; ring
  -- TERM 2: log²(T+2)·x^{σ₀+1}/T² ≤ 12e·X
  have ht2 : (Real.log (Real.exp u + 2)) ^ 2 * x ^ (1 + 1 / Real.log x + 1)
      / (Real.exp u) ^ 2 ≤ 12 * Real.exp 1 * X := by
    have hT2 : (Real.exp u) ^ 2 = Real.exp (2 * u) := by
      rw [sq, ← Real.exp_add]; congr 1; ring
    have hstep : (Real.log (Real.exp u + 2)) ^ 2 * x ^ (1 + 1 / Real.log x + 1)
        / (Real.exp u) ^ 2
        = (Real.log (Real.exp u + 2)) ^ 2 * Real.exp (-(2 * u))
          * (Real.exp 1 * x ^ (2 : ℝ)) := by
      rw [hxσ, hT2, div_eq_mul_inv, ← Real.exp_neg]
      ring
    rw [hstep]
    have hLsq : (Real.log (Real.exp u + 2)) ^ 2 ≤ (u + 2) ^ 2 := by
      nlinarith [hLpos.le, hLhi]
    calc (Real.log (Real.exp u + 2)) ^ 2 * Real.exp (-(2 * u))
          * (Real.exp 1 * x ^ (2 : ℝ))
        ≤ (u + 2) ^ 2 * Real.exp (-(2 * u)) * (Real.exp 1 * x ^ (2 : ℝ)) := by
          apply mul_le_mul_of_nonneg_right _ (by positivity)
          exact mul_le_mul_of_nonneg_right hLsq (Real.exp_pos _).le
      _ ≤ 12 * Real.exp (-u) * (Real.exp 1 * x ^ (2 : ℝ)) :=
          mul_le_mul_of_nonneg_right (sq_shift_absorb hu0) (by positivity)
      _ ≤ 12 * Real.exp (-(c * u)) * (Real.exp 1 * x ^ (2 : ℝ)) := by
          apply mul_le_mul_of_nonneg_right _ (by positivity)
          apply mul_le_mul_of_nonneg_left _ (by norm_num)
          apply Real.exp_le_exp.mpr
          have h1 := mul_le_mul_of_nonneg_right hc_le2 hu0
          linarith
      _ = 12 * Real.exp 1 * X := by rw [hXdef]; ring
  -- TERM 3: ‖L(Λ)(σ₀)‖·x^{σ₀+1}/T ≤ (16+A)e·X
  have ht3 : ‖LSeries (↗Λ) (((1 + 1 / Real.log x : ℝ)) : ℂ)‖
      * x ^ (1 + 1 / Real.log x + 1) / Real.exp u
      ≤ (16 + A) * Real.exp 1 * X := by
    have hnorm := hA (1 + 1 / Real.log x) hσ₀1 hσ₀2
    have hinv : 1 / (1 + 1 / Real.log x - 1) = Real.log x := by
      rw [add_sub_cancel_left, one_div_one_div]
    rw [hinv] at hnorm
    have hnorm2 : ‖LSeries (↗Λ) (((1 + 1 / Real.log x : ℝ)) : ℂ)‖ ≤ u ^ 2 + A := by
      rw [hu2]; exact hnorm
    have hstep : ‖LSeries (↗Λ) (((1 + 1 / Real.log x : ℝ)) : ℂ)‖
        * x ^ (1 + 1 / Real.log x + 1) / Real.exp u
        = ‖LSeries (↗Λ) (((1 + 1 / Real.log x : ℝ)) : ℂ)‖ * Real.exp (-u)
          * (Real.exp 1 * x ^ (2 : ℝ)) := by
      rw [hxσ, div_eq_mul_inv, ← Real.exp_neg]
      ring
    rw [hstep]
    calc ‖LSeries (↗Λ) (((1 + 1 / Real.log x : ℝ)) : ℂ)‖ * Real.exp (-u)
          * (Real.exp 1 * x ^ (2 : ℝ))
        ≤ (u ^ 2 + A) * Real.exp (-u) * (Real.exp 1 * x ^ (2 : ℝ)) := by
          apply mul_le_mul_of_nonneg_right _ (by positivity)
          exact mul_le_mul_of_nonneg_right hnorm2 (Real.exp_pos _).le
      _ ≤ (16 + A) * Real.exp (-(u / 2)) * (Real.exp 1 * x ^ (2 : ℝ)) :=
          mul_le_mul_of_nonneg_right (sq_plus_const_absorb hA0 hu0) (by positivity)
      _ ≤ (16 + A) * Real.exp (-(c * u)) * (Real.exp 1 * x ^ (2 : ℝ)) := by
          apply mul_le_mul_of_nonneg_right _ (by positivity)
          apply mul_le_mul_of_nonneg_left _ (by linarith)
          apply Real.exp_le_exp.mpr
          have h1 := mul_le_mul_of_nonneg_right hc_le2 hu0
          linarith
      _ = (16 + A) * Real.exp 1 * X := by rw [hXdef]; ring
  -- assemble
  calc ‖(2 * Real.pi : ℝ) • psi1C x - 2 * (Real.pi : ℂ) * ((x : ℂ) ^ (2 : ℂ) / 2)‖
      ≤ K * ((1 + Real.log (Real.exp u + 2))
            * x ^ (1 - c₁ / Real.log (Real.exp u + 2) + 1)
          + (Real.log (Real.exp u + 2)) ^ 2 * x ^ (1 + 1 / Real.log x + 1)
            / (Real.exp u) ^ 2
          + ‖LSeries (↗Λ) (((1 + 1 / Real.log x : ℝ)) : ℂ)‖
            * x ^ (1 + 1 / Real.log x + 1) / Real.exp u) := hmain'
    _ ≤ K * ((Real.exp (2 * c₁) * (3 + 2 / c₁) + 12 * Real.exp 1
          + (16 + A) * Real.exp 1) * X) := by
        apply mul_le_mul_of_nonneg_left _ hK.le
        linarith [ht1, ht2, ht3]
    _ = K * (Real.exp (2 * c₁) * (3 + 2 / c₁) + 12 * Real.exp 1
          + (16 + A) * Real.exp 1) * x ^ (2 : ℝ) * Real.exp (-(c * u)) := by
        rw [hXdef]; ring

#print axioms smoothed_pnt_dlvp
