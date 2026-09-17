import LiPartsIdentity
import PNTThetaError

set_option maxHeartbeats 4000000

open Chebyshev Real MeasureTheory intervalIntegral

/-!
# §1706b ★★★ — THE ERROR-INTEGRAL BOUND: `|∫₂ˣ (θt−t)/(t log²t) dt| ≤ (7+3C)·x·e^{-(c/2)√log x}`

Split at `y = √x`.  Below `y`: `|θt − t| ≤ t` (Chebyshev's `θ ≤ log 4·x`), so the integrand is
`≤ 1/log²t ≤ 3` and the chunk is `≤ 3√x ≤ 3x·e^{-(c/2)√log x}` (since `√L ≥ c`).  Above `y`:
the θ-error bound gives integrand `≤ 3C·e^{-c√log y}` and `log y = L/2`, `√(L/2) ≥ √L/2`, so the
chunk is `≤ 3C·x·e^{-(c/2)√log x}`.  Parametrized by the §1705 constants — assembled in §1706.
-/

theorem one_div_log_sq_le_three {t : ℝ} (ht : 2 ≤ t) : 1 / Real.log t ^ 2 ≤ 3 := by
  have h2 : (0.6931471803 : ℝ) < Real.log 2 := Real.log_two_gt_d9
  have hmono : Real.log 2 ≤ Real.log t := Real.log_le_log (by norm_num) ht
  have hpos : 0 < Real.log t := by linarith
  rw [div_le_iff₀ (by positivity)]
  nlinarith

theorem abs_g_le_of_theta_bound {t B : ℝ} (ht : 2 ≤ t) (hB : 0 ≤ B)
    (hθ : |θ t - t| ≤ B * t) :
    |θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2| ≤ 3 * B := by
  have ht0 : (0 : ℝ) < t := by linarith
  have h2 : (0.6931471803 : ℝ) < Real.log 2 := Real.log_two_gt_d9
  have hmono : Real.log 2 ≤ Real.log t := Real.log_le_log (by norm_num) ht
  have hlog : 0 < Real.log t := by linarith
  have hden : 0 < t * Real.log t ^ 2 := by positivity
  have heq : θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2
      = (θ t - t) / (t * Real.log t ^ 2) := by
    field_simp
  rw [heq, abs_div, abs_of_pos hden, div_le_iff₀ hden]
  have h13 : 1 / Real.log t ^ 2 ≤ 3 := one_div_log_sq_le_three ht
  have hlog2ge : 1 ≤ 3 * Real.log t ^ 2 := by
    rw [div_le_iff₀ (by positivity)] at h13
    linarith
  calc |θ t - t| ≤ B * t := hθ
    _ ≤ B * t * (3 * Real.log t ^ 2) := by nlinarith [mul_nonneg hB ht0.le]
    _ = 3 * B * (t * Real.log t ^ 2) := by ring

theorem theta_error_integral_bound {C c x₀θ x : ℝ} (hC : 0 < C) (hc : 0 < c)
    (hx₀θ : 2 ≤ x₀θ)
    (hbound : ∀ t : ℝ, x₀θ ≤ t →
      |θ t - t| ≤ C * t * Real.exp (-(c * Real.sqrt (Real.log t))))
    (hx : x₀θ ^ 2 ≤ x) (hx4 : 4 ≤ x) (hxc : Real.exp (c ^ 2) ≤ x) :
    |∫ t in (2:ℝ)..x, (θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2)|
      ≤ (7 + 3 * C) * x * Real.exp (-(c / 2 * Real.sqrt (Real.log x))) := by
  have hx0 : (0 : ℝ) < x := by linarith
  have hx2 : (2 : ℝ) ≤ x := by linarith
  set y := Real.sqrt x with hydef
  have hy0 : 0 ≤ y := Real.sqrt_nonneg _
  have hy2 : (2 : ℝ) ≤ y := by
    have h1 : ((2 : ℝ)) ^ 2 ≤ x := by nlinarith
    have h2 := Real.sqrt_le_sqrt h1
    rwa [Real.sqrt_sq (by norm_num : (0:ℝ) ≤ 2)] at h2
  have hyx : y ≤ x := by
    have h1 : x ≤ x ^ 2 := by nlinarith
    have h2 := Real.sqrt_le_sqrt h1
    rwa [Real.sqrt_sq hx0.le] at h2
  have hyθ : x₀θ ≤ y := by
    have h2 := Real.sqrt_le_sqrt hx
    rwa [Real.sqrt_sq (by linarith : (0:ℝ) ≤ x₀θ)] at h2
  have hL0 : (0 : ℝ) ≤ Real.log x := Real.log_nonneg (by linarith)
  have hLc : c ^ 2 ≤ Real.log x := (Real.le_log_iff_exp_le hx0).mpr hxc
  have hsqLc : c ≤ Real.sqrt (Real.log x) := by
    have h2 := Real.sqrt_le_sqrt hLc
    rwa [Real.sqrt_sq hc.le] at h2
  -- integrability on [2,x] and the two chunks
  have hIθ : IntervalIntegrable (fun t : ℝ => θ t / (t * Real.log t ^ 2))
      MeasureTheory.volume 2 x := by
    rw [intervalIntegrable_iff, Set.uIoc_of_le hx2]
    exact (integrableOn_theta_div_id_mul_log_sq x).mono_set Set.Ioc_subset_Icc_self
  have hI2 : IntervalIntegrable (fun t : ℝ => 1 / Real.log t ^ 2)
      MeasureTheory.volume 2 x := (continuousOn_one_div_log_sq hx2).intervalIntegrable
  have hIg := hIθ.sub hI2
  have hsub1 : Set.uIcc (2:ℝ) y ⊆ Set.uIcc (2:ℝ) x := by
    rw [Set.uIcc_of_le hy2, Set.uIcc_of_le hx2]
    exact Set.Icc_subset_Icc le_rfl hyx
  have hsub2 : Set.uIcc y x ⊆ Set.uIcc (2:ℝ) x := by
    rw [Set.uIcc_of_le hyx, Set.uIcc_of_le hx2]
    exact Set.Icc_subset_Icc hy2 le_rfl
  have hIg1 := hIg.mono_set hsub1
  have hIg2 := hIg.mono_set hsub2
  -- pointwise bounds
  have hpt1 : ∀ t ∈ Set.Icc (2:ℝ) y,
      |θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2| ≤ 3 := by
    intro t ht
    have ht2 : (2:ℝ) ≤ t := ht.1
    have ht0 : (0:ℝ) < t := by linarith
    have hθb : |θ t - t| ≤ 1 * t := by
      rw [one_mul, abs_le]
      constructor
      · have := theta_nonneg t
        linarith
      · have h4 := theta_le_log4_mul_x ht0.le
        have hlog4 : Real.log 4 ≤ 2 := by
          have h2 : Real.log 4 = 2 * Real.log 2 := by
            rw [show (4:ℝ) = 2 ^ 2 by norm_num, Real.log_pow]
            push_cast
            ring
          have := Real.log_two_lt_d9
          nlinarith
        nlinarith [theta_nonneg t]
    have := abs_g_le_of_theta_bound ht2 (by norm_num : (0:ℝ) ≤ 1) hθb
    linarith
  have hpt2 : ∀ t ∈ Set.Icc y x,
      |θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2|
        ≤ 3 * (C * Real.exp (-(c * Real.sqrt (Real.log y)))) := by
    intro t ht
    have hty : y ≤ t := ht.1
    have ht2 : (2:ℝ) ≤ t := by linarith
    have hθb : |θ t - t|
        ≤ (C * Real.exp (-(c * Real.sqrt (Real.log t)))) * t := by
      have := hbound t (by linarith)
      calc |θ t - t| ≤ C * t * Real.exp (-(c * Real.sqrt (Real.log t))) := this
        _ = (C * Real.exp (-(c * Real.sqrt (Real.log t)))) * t := by ring
    have hCB : (0:ℝ) ≤ C * Real.exp (-(c * Real.sqrt (Real.log t))) := by positivity
    have h1 := abs_g_le_of_theta_bound ht2 hCB hθb
    have hEmono : Real.exp (-(c * Real.sqrt (Real.log t)))
        ≤ Real.exp (-(c * Real.sqrt (Real.log y))) := by
      apply Real.exp_le_exp.mpr
      have hlm : Real.log y ≤ Real.log t := Real.log_le_log (by linarith) hty
      have hsm := Real.sqrt_le_sqrt hlm
      nlinarith [hsm]
    calc |θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2|
        ≤ 3 * (C * Real.exp (-(c * Real.sqrt (Real.log t)))) := h1
      _ ≤ 3 * (C * Real.exp (-(c * Real.sqrt (Real.log y)))) := by
          have := mul_le_mul_of_nonneg_left hEmono hC.le
          linarith
  -- assemble the chunks
  have habs := intervalIntegral.abs_integral_le_integral_abs (μ := MeasureTheory.volume)
    (f := fun t : ℝ => θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2) hx2
  have hsplit : (∫ t in (2:ℝ)..x, |θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2|)
      = (∫ t in (2:ℝ)..y, |θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2|)
        + ∫ t in y..x, |θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2| :=
    (intervalIntegral.integral_add_adjacent_intervals hIg1.abs hIg2.abs).symm
  have hchunk1 : (∫ t in (2:ℝ)..y, |θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2|)
      ≤ 3 * y := by
    have hmono := intervalIntegral.integral_mono_on (μ := MeasureTheory.volume) hy2
      hIg1.abs (intervalIntegrable_const (c := (3:ℝ))) hpt1
    rw [intervalIntegral.integral_const, smul_eq_mul] at hmono
    nlinarith [hmono]
  have hchunk2 : (∫ t in y..x, |θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2|)
      ≤ 3 * C * Real.exp (-(c * Real.sqrt (Real.log y))) * x := by
    have hmono := intervalIntegral.integral_mono_on (μ := MeasureTheory.volume) hyx
      hIg2.abs (intervalIntegrable_const
        (c := 3 * (C * Real.exp (-(c * Real.sqrt (Real.log y)))))) hpt2
    rw [intervalIntegral.integral_const, smul_eq_mul] at hmono
    have hE0 : (0:ℝ) < Real.exp (-(c * Real.sqrt (Real.log y))) := Real.exp_pos _
    nlinarith [hmono, mul_pos hC hE0]
  -- E(y) at scale x: log y = L/2, √(L/2) ≥ √L/2
  have hEy : Real.exp (-(c * Real.sqrt (Real.log y)))
      ≤ Real.exp (-(c / 2 * Real.sqrt (Real.log x))) := by
    apply Real.exp_le_exp.mpr
    have hlogy : Real.log y = Real.log x / 2 := Real.log_sqrt hx0.le
    have hs2 : Real.sqrt (Real.log x) ≤ 2 * Real.sqrt (Real.log x / 2) := by
      have h1 : Real.log x ≤ 2 * Real.log x := by linarith
      have h2 := Real.sqrt_le_sqrt h1
      have h3 : Real.sqrt (2 * Real.log x) = 2 * Real.sqrt (Real.log x / 2) := by
        rw [show 2 * Real.log x = 4 * (Real.log x / 2) by ring,
          show (4 : ℝ) * (Real.log x / 2) = 2 ^ 2 * (Real.log x / 2) by norm_num,
          Real.sqrt_mul (by norm_num) , Real.sqrt_sq (by norm_num : (0:ℝ) ≤ 2)]
      rw [h3] at h2
      exact h2
    rw [hlogy]
    nlinarith [hs2, hc.le]
  -- √x ≤ x·e^{-(c/2)√L}
  have hsx : y ≤ x * Real.exp (-(c / 2 * Real.sqrt (Real.log x))) := by
    have hyexp : Real.exp (Real.log x / 2) = y := by
      rw [hydef, Real.sqrt_eq_rpow, Real.rpow_def_of_pos hx0]
      congr 1
      ring
    have hxexp : Real.exp (Real.log x) = x := Real.exp_log hx0
    have hkey : Real.log x / 2 ≤ Real.log x + -(c / 2 * Real.sqrt (Real.log x)) := by
      have h1 : c / 2 * Real.sqrt (Real.log x)
          ≤ Real.sqrt (Real.log x) / 2 * Real.sqrt (Real.log x) := by
        have := mul_le_mul_of_nonneg_right hsqLc (Real.sqrt_nonneg (Real.log x))
        nlinarith [this]
      have h2 : Real.sqrt (Real.log x) * Real.sqrt (Real.log x) = Real.log x :=
        Real.mul_self_sqrt hL0
      nlinarith [h1, h2]
    calc y = Real.exp (Real.log x / 2) := hyexp.symm
      _ ≤ Real.exp (Real.log x + -(c / 2 * Real.sqrt (Real.log x))) :=
          Real.exp_le_exp.mpr hkey
      _ = x * Real.exp (-(c / 2 * Real.sqrt (Real.log x))) := by
          rw [Real.exp_add, hxexp]
  -- total
  have hE0 : (0:ℝ) ≤ Real.exp (-(c / 2 * Real.sqrt (Real.log x))) := (Real.exp_pos _).le
  calc |∫ t in (2:ℝ)..x, (θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2)|
      ≤ ∫ t in (2:ℝ)..x, |θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2| := habs
    _ = _ + _ := hsplit
    _ ≤ 3 * y + 3 * C * Real.exp (-(c * Real.sqrt (Real.log y))) * x :=
        add_le_add hchunk1 hchunk2
    _ ≤ 3 * (x * Real.exp (-(c / 2 * Real.sqrt (Real.log x))))
        + 3 * C * Real.exp (-(c / 2 * Real.sqrt (Real.log x))) * x := by
        have h1 := mul_le_mul_of_nonneg_left hsx (by norm_num : (0:ℝ) ≤ 3)
        have h2 := mul_le_mul_of_nonneg_right
          (mul_le_mul_of_nonneg_left hEy (by positivity : (0:ℝ) ≤ 3 * C)) hx0.le
        linarith
    _ ≤ (7 + 3 * C) * x * Real.exp (-(c / 2 * Real.sqrt (Real.log x))) := by
        nlinarith [mul_nonneg hx0.le hE0]

#print axioms theta_error_integral_bound
