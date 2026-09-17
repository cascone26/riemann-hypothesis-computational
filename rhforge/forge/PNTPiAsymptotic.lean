import PNTPiLiError

set_option maxHeartbeats 4000000

open Chebyshev Real MeasureTheory intervalIntegral
open scoped Nat.Prime

/-!
# §1708 ★★★★ — THE TEXTBOOK FORM: `π(x) = x/log x + O(x/log²x)`

Corollary of §1706: `Li(x) − x/log x = ∫₂ˣdt/log²t − 2/log 2` (parts identity §1706a), the
integral is `≤ 7x/log²x` (√x-split: `≤ 3` below, `≤ 4/log²x` above, and `√x·log²x ≤ x` once
`log x ≥ 64` since `e^{L/2} = (e^{L/8})⁴ ≥ (L/8)⁴ ≥ L²`), and the dlVP error is sub-polynomial:
`L²·e^{-c√L} = (L·e^{-(c/2)√L})² ≤ (16/c²)²` by squaring §1703a's absorption at `u = √L`.
-/

theorem integral_inv_log_sq_bounds {x : ℝ} (hx : Real.exp 64 ≤ x) :
    (0 ≤ ∫ t in (2:ℝ)..x, 1 / Real.log t ^ 2)
    ∧ (∫ t in (2:ℝ)..x, 1 / Real.log t ^ 2) ≤ 7 * x / Real.log x ^ 2 := by
  have hx0 : (0 : ℝ) < x := lt_of_lt_of_le (Real.exp_pos _) hx
  have hL64 : (64 : ℝ) ≤ Real.log x := (Real.le_log_iff_exp_le hx0).mpr hx
  have hL0 : (0 : ℝ) < Real.log x := by linarith
  have hx4 : (4 : ℝ) ≤ x := by
    have h1 : Real.exp 64 ≥ 64 + 1 := by have := Real.add_one_le_exp (64:ℝ); linarith
    linarith
  have hx2 : (2 : ℝ) ≤ x := by linarith
  set y := Real.sqrt x with hydef
  have hy2 : (2 : ℝ) ≤ y := by
    have h1 : ((2 : ℝ)) ^ 2 ≤ x := by nlinarith [hx4]
    have h2 := Real.sqrt_le_sqrt h1
    rwa [Real.sqrt_sq (by norm_num : (0:ℝ) ≤ 2)] at h2
  have hyx : y ≤ x := by
    have h1 : x ≤ x ^ 2 := by nlinarith
    have h2 := Real.sqrt_le_sqrt h1
    rwa [Real.sqrt_sq hx0.le] at h2
  have hlogy : Real.log y = Real.log x / 2 := Real.log_sqrt hx0.le
  -- integrability on the pieces
  have hI : IntervalIntegrable (fun t : ℝ => 1 / Real.log t ^ 2)
      MeasureTheory.volume 2 x := (continuousOn_one_div_log_sq hx2).intervalIntegrable
  have hsub1 : Set.uIcc (2:ℝ) y ⊆ Set.uIcc (2:ℝ) x := by
    rw [Set.uIcc_of_le hy2, Set.uIcc_of_le hx2]
    exact Set.Icc_subset_Icc le_rfl hyx
  have hsub2 : Set.uIcc y x ⊆ Set.uIcc (2:ℝ) x := by
    rw [Set.uIcc_of_le hyx, Set.uIcc_of_le hx2]
    exact Set.Icc_subset_Icc hy2 le_rfl
  have hI1 := hI.mono_set hsub1
  have hI2 := hI.mono_set hsub2
  constructor
  · apply intervalIntegral.integral_nonneg hx2
    intro t ht
    positivity
  · have hsplit : (∫ t in (2:ℝ)..x, 1 / Real.log t ^ 2)
        = (∫ t in (2:ℝ)..y, 1 / Real.log t ^ 2) + ∫ t in y..x, 1 / Real.log t ^ 2 :=
      (intervalIntegral.integral_add_adjacent_intervals hI1 hI2).symm
    have hchunk1 : (∫ t in (2:ℝ)..y, 1 / Real.log t ^ 2) ≤ 3 * y := by
      have hmono := intervalIntegral.integral_mono_on (μ := MeasureTheory.volume) hy2
        hI1 (intervalIntegrable_const (c := (3:ℝ)))
        (fun t ht => one_div_log_sq_le_three ht.1)
      rw [intervalIntegral.integral_const, smul_eq_mul] at hmono
      nlinarith [hmono]
    have hchunk2 : (∫ t in y..x, 1 / Real.log t ^ 2)
        ≤ 4 / Real.log x ^ 2 * x := by
      have hpt : ∀ t ∈ Set.Icc y x, 1 / Real.log t ^ 2 ≤ 4 / Real.log x ^ 2 := by
        intro t ht
        have hlt : Real.log y ≤ Real.log t := Real.log_le_log (by linarith) ht.1
        have hly0 : 0 < Real.log y := by rw [hlogy]; linarith
        have hlt0 : 0 < Real.log t := by linarith
        rw [div_le_div_iff₀ (by positivity) (by positivity)]
        have h1 : Real.log x = 2 * Real.log y := by rw [hlogy]; ring
        nlinarith [hlt, hly0]
      have hmono := intervalIntegral.integral_mono_on (μ := MeasureTheory.volume) hyx
        hI2 (intervalIntegrable_const (c := 4 / Real.log x ^ 2)) hpt
      rw [intervalIntegral.integral_const, smul_eq_mul] at hmono
      have h40 : (0:ℝ) ≤ 4 / Real.log x ^ 2 := by positivity
      nlinarith [hmono, mul_nonneg h40 (Real.sqrt_nonneg x)]
    -- 3√x ≤ 3x/log²x  (e^{L/2} ≥ (L/8)⁴ ≥ L²)
    have hsx : y * Real.log x ^ 2 ≤ x := by
      have hyexp : Real.exp (Real.log x / 2) = y := by
        rw [hydef, Real.sqrt_eq_rpow, Real.rpow_def_of_pos hx0]
        congr 1
        ring
      have h8 : 1 + Real.log x / 8 ≤ Real.exp (Real.log x / 8) := by
        have := Real.add_one_le_exp (Real.log x / 8)
        linarith
      have h4 : Real.exp (Real.log x / 2) = Real.exp (Real.log x / 8) ^ 4 := by
        rw [← Real.exp_nat_mul]
        congr 1
        push_cast
        ring
      have hLsq : Real.log x ^ 2 ≤ Real.exp (Real.log x / 2) := by
        rw [h4]
        have hstep : (Real.log x / 8) ^ 4 ≤ Real.exp (Real.log x / 8) ^ 4 := by
          have h0 : (0:ℝ) ≤ Real.log x / 8 := by linarith
          gcongr
          linarith
        have hLL : (4096:ℝ) ≤ Real.log x ^ 2 := by nlinarith [hL64]
        have hpoly : Real.log x ^ 2 ≤ (Real.log x / 8) ^ 4 := by
          nlinarith [mul_le_mul_of_nonneg_right hLL (sq_nonneg (Real.log x))]
        linarith
      calc y * Real.log x ^ 2 ≤ y * Real.exp (Real.log x / 2) := by
            have := mul_le_mul_of_nonneg_left hLsq (Real.sqrt_nonneg x)
            simpa [hydef] using this
        _ = Real.exp (Real.log x / 2) * Real.exp (Real.log x / 2) := by rw [hyexp]
        _ = Real.exp (Real.log x) := by rw [← Real.exp_add]; congr 1; ring
        _ = x := Real.exp_log hx0
    rw [hsplit]
    have h3y : 3 * y ≤ 3 * x / Real.log x ^ 2 := by
      rw [le_div_iff₀ (by positivity)]
      nlinarith [hsx]
    have h4x : 4 / Real.log x ^ 2 * x = 4 * x / Real.log x ^ 2 := by ring
    rw [h4x] at hchunk2
    have : 3 * x / Real.log x ^ 2 + 4 * x / Real.log x ^ 2 = 7 * x / Real.log x ^ 2 := by
      ring
    linarith

theorem sq_log_mul_exp_le {c L : ℝ} (hc : 0 < c) (hL : 0 ≤ L) :
    L ^ 2 * Real.exp (-(c * Real.sqrt L)) ≤ 256 / c ^ 4 := by
  have hu := sq_mul_exp_neg_le (ε := c / 2) (u := Real.sqrt L) (by positivity)
    (Real.sqrt_nonneg L)
  rw [Real.sq_sqrt hL] at hu
  -- hu : L·e^{-(c/2)√L} ≤ 4/(c/2)² = 16/c²
  have hu' : L * Real.exp (-(c / 2 * Real.sqrt L)) ≤ 16 / c ^ 2 := by
    have h1 : (4:ℝ) / (c / 2) ^ 2 = 16 / c ^ 2 := by
      field_simp
      ring
    linarith [hu, h1.symm.le, h1.le]
  have hsq := mul_self_le_mul_self
    (by positivity : (0:ℝ) ≤ L * Real.exp (-(c / 2 * Real.sqrt L))) hu'
  have hexp : Real.exp (-(c / 2 * Real.sqrt L)) * Real.exp (-(c / 2 * Real.sqrt L))
      = Real.exp (-(c * Real.sqrt L)) := by
    rw [← Real.exp_add]
    congr 1
    ring
  have h16 : (16:ℝ) / c ^ 2 * (16 / c ^ 2) = 256 / c ^ 4 := by
    field_simp
    ring
  nlinarith [hsq, hexp, h16]

theorem pnt_pi_asymptotic :
    ∃ C x₀ : ℝ, 0 < C ∧ 2 ≤ x₀ ∧ ∀ x : ℝ, x₀ ≤ x →
      |(π ⌊x⌋₊ : ℝ) - x / Real.log x| ≤ C * x / Real.log x ^ 2 := by
  obtain ⟨C, c, x₀π, hC, hc, hx₀π, hπ⟩ := pnt_pi_li_error
  refine ⟨8 + 256 * C / c ^ 4, max x₀π (Real.exp 64), by positivity,
    le_trans hx₀π (le_max_left _ _), ?_⟩
  intro x hx
  have hxπ : x₀π ≤ x := le_trans (le_max_left _ _) hx
  have hxe : Real.exp 64 ≤ x := le_trans (le_max_right _ _) hx
  have hx0 : (0 : ℝ) < x := lt_of_lt_of_le (Real.exp_pos _) hxe
  have hx2 : (2 : ℝ) ≤ x := by linarith [hx₀π]
  have hL64 : (64 : ℝ) ≤ Real.log x := (Real.le_log_iff_exp_le hx0).mpr hxe
  have hL0 : (0 : ℝ) < Real.log x := by linarith
  obtain ⟨hint0, hintle⟩ := integral_inv_log_sq_bounds hxe
  have hparts := li_parts_identity hx2
  have hπx := hπ x hxπ
  -- |Li − x/log x| ≤ ∫1/log² + 2/log 2 ≤ 7x/L² + x/L²
  have hlog2 : (0.6931471803 : ℝ) < Real.log 2 := Real.log_two_gt_d9
  have hconst : 2 / Real.log 2 ≤ x / Real.log x ^ 2 := by
    have h1 : 2 / Real.log 2 ≤ 3 := by
      rw [div_le_iff₀ (by linarith)]
      linarith
    have h2 : 3 * Real.log x ^ 2 ≤ x := by
      -- e^L ≥ (L/8)⁴·... :  L⁴/4096 ≥ 3L² for L ≥ 111; L ≥ 64: (L/8)⁴ = L⁴/4096 ≥ 3L² ⟺ L² ≥ 12288, ok for L ≥ 111... use e^L = (e^{L/2})² ≥ (L²)² since e^{L/2} ≥ L² at L ≥ 64
      have h8 : 1 + Real.log x / 8 ≤ Real.exp (Real.log x / 8) := by
        have := Real.add_one_le_exp (Real.log x / 8)
        linarith
      have h4 : Real.exp (Real.log x / 2) = Real.exp (Real.log x / 8) ^ 4 := by
        rw [← Real.exp_nat_mul]
        congr 1
        push_cast
        ring
      have hLsq : Real.log x ^ 2 ≤ Real.exp (Real.log x / 2) := by
        rw [h4]
        have hstep : (Real.log x / 8) ^ 4 ≤ Real.exp (Real.log x / 8) ^ 4 := by
          have h0 : (0:ℝ) ≤ Real.log x / 8 := by linarith
          gcongr
          linarith
        have hLL : (4096:ℝ) ≤ Real.log x ^ 2 := by nlinarith [hL64]
        have hpoly : Real.log x ^ 2 ≤ (Real.log x / 8) ^ 4 := by
          nlinarith [mul_le_mul_of_nonneg_right hLL (sq_nonneg (Real.log x))]
        linarith
      have hxeq : x = Real.exp (Real.log x / 2) * Real.exp (Real.log x / 2) := by
        rw [← Real.exp_add, show Real.log x / 2 + Real.log x / 2 = Real.log x by ring,
          Real.exp_log hx0]
      have hE2 : Real.exp (Real.log x / 2) ≥ 3 := by
        have : (3:ℝ) ≤ 1 + Real.log x / 8 := by linarith
        have h4' : Real.exp (Real.log x / 8) ≤ Real.exp (Real.log x / 2) :=
          Real.exp_le_exp.mpr (by linarith)
        linarith [h8]
      nlinarith [hLsq, hE2, Real.exp_pos (Real.log x / 2)]
    rw [div_le_div_iff₀ (by linarith) (by positivity)]
    nlinarith [h2]
  -- dlVP error is sub-polynomial: x·E = x/L²·(L²E) ≤ (256/c⁴)·x/L²
  have hsub : C * x * Real.exp (-(c * Real.sqrt (Real.log x)))
      ≤ 256 * C / c ^ 4 * x / Real.log x ^ 2 := by
    have h1 := sq_log_mul_exp_le (L := Real.log x) hc hL0.le
    have hE0 : (0:ℝ) < Real.exp (-(c * Real.sqrt (Real.log x))) := Real.exp_pos _
    calc C * x * Real.exp (-(c * Real.sqrt (Real.log x)))
        = C * x * (Real.log x ^ 2 * Real.exp (-(c * Real.sqrt (Real.log x))))
          / Real.log x ^ 2 := by
          field_simp
      _ ≤ C * x * (256 / c ^ 4) / Real.log x ^ 2 := by
          gcongr
      _ = 256 * C / c ^ 4 * x / Real.log x ^ 2 := by ring
  -- assemble:  π − x/log x = (π − Li) + (∫1/log² − 2/log 2)
  have hdecomp : (π ⌊x⌋₊ : ℝ) - x / Real.log x
      = ((π ⌊x⌋₊ : ℝ) - ∫ t in (2:ℝ)..x, 1 / Real.log t)
        + ((∫ t in (2:ℝ)..x, 1 / Real.log t ^ 2) - 2 / Real.log 2) := by
    linarith [hparts]
  rw [hdecomp]
  have habsB : |(∫ t in (2:ℝ)..x, 1 / Real.log t ^ 2) - 2 / Real.log 2|
      ≤ (∫ t in (2:ℝ)..x, 1 / Real.log t ^ 2) + 2 / Real.log 2 := by
    rw [abs_le]
    constructor
    · have h20 : (0:ℝ) < 2 / Real.log 2 := by positivity
      linarith
    · have h20 : (0:ℝ) < 2 / Real.log 2 := by positivity
      linarith
  calc |((π ⌊x⌋₊ : ℝ) - ∫ t in (2:ℝ)..x, 1 / Real.log t)
        + ((∫ t in (2:ℝ)..x, 1 / Real.log t ^ 2) - 2 / Real.log 2)|
      ≤ |(π ⌊x⌋₊ : ℝ) - ∫ t in (2:ℝ)..x, 1 / Real.log t|
        + |(∫ t in (2:ℝ)..x, 1 / Real.log t ^ 2) - 2 / Real.log 2| := abs_add_le _ _
    _ ≤ C * x * Real.exp (-(c * Real.sqrt (Real.log x)))
        + ((∫ t in (2:ℝ)..x, 1 / Real.log t ^ 2) + 2 / Real.log 2) :=
        add_le_add hπx habsB
    _ ≤ 256 * C / c ^ 4 * x / Real.log x ^ 2
        + (7 * x / Real.log x ^ 2 + x / Real.log x ^ 2) := by
        linarith [hsub, hintle, hconst]
    _ = (8 + 256 * C / c ^ 4) * x / Real.log x ^ 2 := by ring

#print axioms pnt_pi_asymptotic
