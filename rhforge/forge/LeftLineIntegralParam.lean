import LeftLinePackage
import LeftLineIntegrable2
import LogSqTailIntegral

set_option maxHeartbeats 4000000

open Complex Metric intervalIntegral

/-!
# §1697a — the §1696 bound, PARAMETRIZED by the package data (for coherent composition): `‖∫ F‖ ≤ K·(1+log(T+2))·x^{σ₁(T)+1}`

The Perron error term, quantified: the full truncated integral of `F` along the zero-free line
`Re = σ₁(T) = 1 - c₁/log(T+2)` is at most `K(1+log(T+2))·x^{σ₁(T)+1}` — uniformly in `T`. Three
pieces: outer `|t| ∈ [Tsp, T]` via the `t`-local `A·log²(|t|+2)/t²` majorant and the §1693
convergent tail (`≤ 96A` each side, the left one by `t ↦ -t`); middle `|t| ≤ Tsp` via the pole-split
constant `+ L/c₁` over a fixed-length interval. Since `x^{σ₁(T)+1} = x²·x^{-c₁/log(T+2)}`, choosing
`log T ≍ √log x` later gives the classical `x²·e^{-c√log x}` error of smoothed PNT.
-/

theorem left_line_integral_bound_param {c₁ A T₁ Tsp C : ℝ}
    (hc₁ : 0 < c₁) (hA : 0 < A) (hT₁6 : 6 ≤ T₁) (hTsp6 : 6 ≤ Tsp) (hTspT₁ : Tsp ≤ T₁)
    (hC : 0 ≤ C)
    (hpkg : ∀ T : ℝ, T₁ ≤ T →
      (∀ σ : ℝ, 1 - c₁ / Real.log (T + 2) ≤ σ → σ < 1 → ∀ t : ℝ, |t| ≤ T →
        (Tsp ≤ |t| →
          ‖logDeriv riemannZeta (((σ : ℝ) : ℂ) + (t : ℝ) * Complex.I)‖
            ≤ A * (Real.log (|t| + 2)) ^ 2) ∧
        (|t| ≤ Tsp →
          ‖logDeriv riemannZeta (((σ : ℝ) : ℂ) + (t : ℝ) * Complex.I)‖
            ≤ C + 1 / (1 - σ))) ∧
      (∀ z : ℂ, 1 - c₁ / Real.log (T + 2) ≤ z.re → z.re < 1 → |z.im| ≤ T + 1 →
        riemannZeta z ≠ 0)) :
    ∀ T : ℝ, max T₁ (Real.exp (2 * c₁)) ≤ T → ∀ x : ℝ, 1 ≤ x →
      ‖∫ t in (-T)..T,
          -logDeriv riemannZeta
              (((1 - c₁ / Real.log (T + 2) : ℝ) : ℂ) + (t : ℝ) * Complex.I) *
              (x : ℂ) ^ ((((1 - c₁ / Real.log (T + 2) : ℝ) : ℂ) + (t : ℝ) * Complex.I) + 1) /
            ((((1 - c₁ / Real.log (T + 2) : ℝ) : ℂ) + (t : ℝ) * Complex.I) *
              ((((1 - c₁ / Real.log (T + 2) : ℝ) : ℂ) + (t : ℝ) * Complex.I) + 1))‖
        ≤ (192 * A + 4 * Tsp * C + 4 * Tsp / c₁) * (1 + Real.log (T + 2))
          * x ^ (1 - c₁ / Real.log (T + 2) + 1) := by
  intro T hT x hx
  have hTT₁ : T₁ ≤ T := le_trans (le_max_left _ _) hT
  have hT6 : (6 : ℝ) ≤ T := le_trans hT₁6 hTT₁
  have hL : 0 < Real.log (T + 2) := Real.log_pos (by linarith)
  have hL2c : 2 * c₁ ≤ Real.log (T + 2) := by
    rw [← Real.log_exp (2 * c₁)]
    apply Real.log_le_log (Real.exp_pos _)
    have := le_trans (le_max_right T₁ (Real.exp (2 * c₁))) hT
    linarith
  set L : ℝ := Real.log (T + 2) with hLdef
  set σ₁ : ℝ := 1 - c₁ / L with hσ₁def
  have hσ₁half : 1 / 2 ≤ σ₁ := by
    rw [hσ₁def]
    have : c₁ / L ≤ 1 / 2 := by
      rw [div_le_iff₀ hL]
      linarith
    linarith
  have hσ₁lt : σ₁ < 1 := by
    rw [hσ₁def]
    have : 0 < c₁ / L := by positivity
    linarith
  obtain ⟨hbounds, hzf⟩ := hpkg T hTT₁
  -- shorthand
  set s : ℝ → ℂ := fun t => ((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I with hsdef
  set F : ℝ → ℂ := fun t =>
    -logDeriv riemannZeta (s t) * (x : ℂ) ^ (s t + 1) / (s t * (s t + 1)) with hFdef
  have hsre : ∀ t, (s t).re = σ₁ := by intro t; rw [hsdef]; simp
  have hsim : ∀ t, (s t).im = t := by intro t; rw [hsdef]; simp
  -- the x-power factor
  have hxppos : (0 : ℝ) < x ^ (σ₁ + 1) := Real.rpow_pos_of_pos (by linarith) _
  have hxp : ∀ t : ℝ, ‖(x : ℂ) ^ (s t + 1)‖ = x ^ (σ₁ + 1) := by
    intro t
    rw [Complex.norm_cpow_eq_rpow_re_of_pos (by linarith)]
    congr 1
    rw [Complex.add_re, hsre, Complex.one_re]
  -- norm decomposition
  have hFnorm : ∀ t : ℝ, ‖F t‖
      = ‖logDeriv riemannZeta (s t)‖ * x ^ (σ₁ + 1) / ‖s t * (s t + 1)‖ := by
    intro t
    rw [hFdef]
    simp only []
    rw [norm_div, norm_mul, norm_neg, hxp]
  -- denominator bounds
  have hden_large : ∀ t : ℝ, Tsp ≤ |t| → t ^ 2 ≤ ‖s t * (s t + 1)‖ := by
    intro t htb
    have htpos : 0 < |t| := lt_of_lt_of_le (by linarith) htb
    rw [norm_mul]
    have h1 : |t| ≤ ‖s t‖ := by
      calc |t| = |(s t).im| := by rw [hsim]
        _ ≤ ‖s t‖ := Complex.abs_im_le_norm _
    have h2 : |t| ≤ ‖s t + 1‖ := by
      have him1 : (s t + 1).im = t := by rw [Complex.add_im, hsim, Complex.one_im, add_zero]
      calc |t| = |(s t + 1).im| := by rw [him1]
        _ ≤ ‖s t + 1‖ := Complex.abs_im_le_norm _
    calc t ^ 2 = |t| * |t| := by rw [← sq_abs]; ring
      _ ≤ ‖s t‖ * ‖s t + 1‖ := mul_le_mul h1 h2 htpos.le (norm_nonneg _)
  have hden_small : ∀ t : ℝ, 1 / 2 ≤ ‖s t * (s t + 1)‖ := by
    intro t
    rw [norm_mul]
    have h1 : σ₁ ≤ ‖s t‖ := by
      calc σ₁ = |σ₁| := (abs_of_pos (by linarith)).symm
        _ = |(s t).re| := by rw [hsre]
        _ ≤ ‖s t‖ := Complex.abs_re_le_norm _
    have h2 : (1 : ℝ) ≤ ‖s t + 1‖ := by
      have hre : (s t + 1).re = σ₁ + 1 := by rw [Complex.add_re, hsre, Complex.one_re]
      calc (1 : ℝ) ≤ σ₁ + 1 := by linarith
        _ = |σ₁ + 1| := (abs_of_pos (by linarith)).symm
        _ = |(s t + 1).re| := by rw [hre]
        _ ≤ ‖s t + 1‖ := Complex.abs_re_le_norm _
    calc (1 : ℝ) / 2 = (1 / 2) * 1 := by ring
      _ ≤ ‖s t‖ * ‖s t + 1‖ := by
          apply mul_le_mul (le_trans (by linarith) h1) h2 (by norm_num) (norm_nonneg _)
  -- integrability of F on subintervals of [-T, T]
  have hInt : ∀ a b : ℝ, -T ≤ a → a ≤ b → b ≤ T →
      IntervalIntegrable F MeasureTheory.volume a b := by
    intro a b hga hab hbT
    rw [hFdef, hsdef]
    apply left_line_F_intervalIntegrable' hc₁ (by linarith : (2:ℝ) ≤ T) hx
    · rw [← hLdef]
      have h := hσ₁half
      rw [hσ₁def] at h
      linarith
    · intro z hzre hzim
      apply hzf z (le_of_eq hzre.symm) _ hzim
      rw [hzre, ← hLdef, ← hσ₁def]
      exact hσ₁lt
    · intro t ht
      rw [Set.uIcc_of_le hab, Set.mem_Icc] at ht
      rw [abs_le]
      constructor <;> linarith [ht.1, ht.2]
  have hIntL := hInt (-T) (-Tsp) (le_refl _) (by linarith) (by linarith)
  have hIntM := hInt (-Tsp) Tsp (by linarith) (by linarith) (by linarith)
  have hIntR := hInt Tsp T (by linarith) (by linarith) (le_refl _)
  have hIntLM := hIntL.trans hIntM
  -- split the integral
  have hsplit : (∫ t in (-T)..T, F t)
      = ((∫ t in (-T)..(-Tsp), F t) + ∫ t in (-Tsp)..Tsp, F t) + ∫ t in Tsp..T, F t := by
    rw [integral_add_adjacent_intervals hIntL hIntM,
      integral_add_adjacent_intervals hIntLM hIntR]
  -- pointwise majorants
  have hpt_right : ∀ t ∈ Set.Icc Tsp T,
      ‖F t‖ ≤ A * x ^ (σ₁ + 1) * ((Real.log (t + 2)) ^ 2 / t ^ 2) := by
    intro t ht
    have htpos : 0 < t := by linarith [ht.1]
    have habs : |t| = t := abs_of_pos htpos
    have hld := (hbounds σ₁ (by rw [hσ₁def, ← hLdef]) hσ₁lt t
      (by rw [habs]; exact ht.2)).1 (by rw [habs]; exact ht.1)
    rw [habs] at hld
    have hdl := hden_large t (by rw [habs]; exact ht.1)
    have hdenpos : (0 : ℝ) < ‖s t * (s t + 1)‖ :=
      lt_of_lt_of_le (by positivity) hdl
    rw [hFnorm]
    have step1 : ‖logDeriv riemannZeta (s t)‖ * x ^ (σ₁ + 1) / ‖s t * (s t + 1)‖
        ≤ A * (Real.log (t + 2)) ^ 2 * x ^ (σ₁ + 1) / ‖s t * (s t + 1)‖ := by
      gcongr
    have step2 : A * (Real.log (t + 2)) ^ 2 * x ^ (σ₁ + 1) / ‖s t * (s t + 1)‖
        ≤ A * (Real.log (t + 2)) ^ 2 * x ^ (σ₁ + 1) / t ^ 2 :=
      div_le_div_of_nonneg_left
        (mul_nonneg (mul_nonneg hA.le (sq_nonneg _)) hxppos.le)
        (pow_pos htpos 2) hdl
    calc ‖logDeriv riemannZeta (s t)‖ * x ^ (σ₁ + 1) / ‖s t * (s t + 1)‖
        ≤ A * (Real.log (t + 2)) ^ 2 * x ^ (σ₁ + 1) / t ^ 2 := le_trans step1 step2
      _ = A * x ^ (σ₁ + 1) * ((Real.log (t + 2)) ^ 2 / t ^ 2) := by ring
  have hpt_left : ∀ t ∈ Set.Icc (-T) (-Tsp),
      ‖F t‖ ≤ A * x ^ (σ₁ + 1) * ((Real.log (-t + 2)) ^ 2 / t ^ 2) := by
    intro t ht
    have htneg : t < 0 := by linarith [ht.2]
    have habs : |t| = -t := abs_of_neg htneg
    have hld := (hbounds σ₁ (by rw [hσ₁def, ← hLdef]) hσ₁lt t
      (by rw [habs]; linarith [ht.1])).1 (by rw [habs]; linarith [ht.2])
    rw [habs] at hld
    have hdl := hden_large t (by rw [habs]; linarith [ht.2])
    rw [hFnorm]
    have step1 : ‖logDeriv riemannZeta (s t)‖ * x ^ (σ₁ + 1) / ‖s t * (s t + 1)‖
        ≤ A * (Real.log (-t + 2)) ^ 2 * x ^ (σ₁ + 1) / ‖s t * (s t + 1)‖ := by
      gcongr
    have step2 : A * (Real.log (-t + 2)) ^ 2 * x ^ (σ₁ + 1) / ‖s t * (s t + 1)‖
        ≤ A * (Real.log (-t + 2)) ^ 2 * x ^ (σ₁ + 1) / t ^ 2 :=
      div_le_div_of_nonneg_left
        (mul_nonneg (mul_nonneg hA.le (sq_nonneg _)) hxppos.le)
        (pow_two_pos_of_ne_zero (ne_of_lt htneg)) hdl
    calc ‖logDeriv riemannZeta (s t)‖ * x ^ (σ₁ + 1) / ‖s t * (s t + 1)‖
        ≤ A * (Real.log (-t + 2)) ^ 2 * x ^ (σ₁ + 1) / t ^ 2 := le_trans step1 step2
      _ = A * x ^ (σ₁ + 1) * ((Real.log (-t + 2)) ^ 2 / t ^ 2) := by ring
  have hpt_mid : ∀ t ∈ Set.uIoc (-Tsp) Tsp,
      ‖F t‖ ≤ (C + L / c₁) * x ^ (σ₁ + 1) * 2 := by
    intro t ht
    have ht' : t ∈ Set.Icc (-Tsp) Tsp := by
      rw [Set.uIoc_of_le (by linarith)] at ht
      exact ⟨ht.1.le, ht.2⟩
    have habs : |t| ≤ Tsp := by
      rw [abs_le]
      exact ⟨ht'.1, ht'.2⟩
    have hld := (hbounds σ₁ (by rw [hσ₁def, ← hLdef]) hσ₁lt t
      (by linarith [habs, hTspT₁, hTT₁])).2 habs
    have hpole : 1 / (1 - σ₁) = L / c₁ := by
      rw [hσ₁def, show (1 : ℝ) - (1 - c₁ / L) = c₁ / L by ring, one_div_div]
    rw [hpole] at hld
    have hds := hden_small t
    rw [hFnorm]
    have hdenpos : (0 : ℝ) < ‖s t * (s t + 1)‖ := lt_of_lt_of_le (by norm_num) hds
    have step1 : ‖logDeriv riemannZeta (s t)‖ * x ^ (σ₁ + 1) / ‖s t * (s t + 1)‖
        ≤ (C + L / c₁) * x ^ (σ₁ + 1) / ‖s t * (s t + 1)‖ := by
      gcongr
    have step2 : (C + L / c₁) * x ^ (σ₁ + 1) / ‖s t * (s t + 1)‖
        ≤ (C + L / c₁) * x ^ (σ₁ + 1) / (1 / 2) :=
      div_le_div_of_nonneg_left
        (mul_nonneg (add_nonneg hC (by positivity)) hxppos.le)
        (by norm_num) hds
    calc ‖logDeriv riemannZeta (s t)‖ * x ^ (σ₁ + 1) / ‖s t * (s t + 1)‖
        ≤ (C + L / c₁) * x ^ (σ₁ + 1) / (1 / 2) := le_trans step1 step2
      _ = (C + L / c₁) * x ^ (σ₁ + 1) * 2 := by ring
  -- majorant integrability (outer pieces)
  have hmajR : IntervalIntegrable
      (fun t : ℝ => A * x ^ (σ₁ + 1) * ((Real.log (t + 2)) ^ 2 / t ^ 2))
      MeasureTheory.volume Tsp T := by
    apply ContinuousOn.intervalIntegrable
    rw [Set.uIcc_of_le (by linarith)]
    apply ContinuousOn.mul continuousOn_const
    apply ContinuousOn.div
    · apply ContinuousOn.pow
      apply ContinuousOn.log (by fun_prop)
      intro u hu
      have := hu.1
      intro h
      linarith
    · fun_prop
    · intro u hu
      have h1 : Tsp ≤ u := hu.1
      have h0 : (0 : ℝ) < u := by linarith
      positivity
  have hmajL : IntervalIntegrable
      (fun t : ℝ => A * x ^ (σ₁ + 1) * ((Real.log (-t + 2)) ^ 2 / t ^ 2))
      MeasureTheory.volume (-T) (-Tsp) := by
    apply ContinuousOn.intervalIntegrable
    rw [Set.uIcc_of_le (by linarith)]
    apply ContinuousOn.mul continuousOn_const
    apply ContinuousOn.div
    · apply ContinuousOn.pow
      apply ContinuousOn.log (by fun_prop)
      intro u hu
      have := hu.2
      intro h
      linarith
    · fun_prop
    · intro u hu
      have h1 : u ≤ -Tsp := hu.2
      have h0 : u < 0 := by linarith
      have : u ≠ 0 := ne_of_lt h0
      positivity
  -- piece bounds
  have hI3 : ‖∫ t in Tsp..T, F t‖ ≤ 96 * A * x ^ (σ₁ + 1) := by
    calc ‖∫ t in Tsp..T, F t‖
        ≤ ∫ t in Tsp..T, ‖F t‖ :=
          intervalIntegral.norm_integral_le_integral_norm (by linarith)
      _ ≤ ∫ t in Tsp..T, A * x ^ (σ₁ + 1) * ((Real.log (t + 2)) ^ 2 / t ^ 2) :=
          intervalIntegral.integral_mono_on (by linarith) hIntR.norm hmajR hpt_right
      _ = A * x ^ (σ₁ + 1) * ∫ t in Tsp..T, (Real.log (t + 2)) ^ 2 / t ^ 2 :=
          intervalIntegral.integral_const_mul _ _
      _ ≤ A * x ^ (σ₁ + 1) * (96 / Real.sqrt Tsp) := by
          gcongr
          exact integral_log_sq_div_sq_le (by linarith) (by linarith)
      _ ≤ A * x ^ (σ₁ + 1) * 96 := by
          have hs1 : (1 : ℝ) ≤ Real.sqrt Tsp := by
            rw [show (1 : ℝ) = Real.sqrt 1 by rw [Real.sqrt_one]]
            exact Real.sqrt_le_sqrt (by linarith)
          gcongr
          rw [div_le_iff₀ (by linarith)]
          nlinarith
      _ = 96 * A * x ^ (σ₁ + 1) := by ring
  have hI1 : ‖∫ t in (-T)..(-Tsp), F t‖ ≤ 96 * A * x ^ (σ₁ + 1) := by
    have hcomp : (∫ t in (-T)..(-Tsp), (Real.log (-t + 2)) ^ 2 / t ^ 2)
        = ∫ u in Tsp..T, (Real.log (u + 2)) ^ 2 / u ^ 2 := by
      have := intervalIntegral.integral_comp_neg
        (fun u : ℝ => (Real.log (u + 2)) ^ 2 / u ^ 2) (a := -T) (b := -Tsp)
      simp only [neg_neg] at this
      calc (∫ t in (-T)..(-Tsp), (Real.log (-t + 2)) ^ 2 / t ^ 2)
          = ∫ t in (-T)..(-Tsp),
              (fun u : ℝ => (Real.log (u + 2)) ^ 2 / u ^ 2) (-t) := by
            apply intervalIntegral.integral_congr
            intro t _
            simp only []
            rw [neg_pow]
            norm_num
        _ = ∫ u in Tsp..T, (Real.log (u + 2)) ^ 2 / u ^ 2 := this
    calc ‖∫ t in (-T)..(-Tsp), F t‖
        ≤ ∫ t in (-T)..(-Tsp), ‖F t‖ :=
          intervalIntegral.norm_integral_le_integral_norm (by linarith)
      _ ≤ ∫ t in (-T)..(-Tsp), A * x ^ (σ₁ + 1) * ((Real.log (-t + 2)) ^ 2 / t ^ 2) :=
          intervalIntegral.integral_mono_on (by linarith) hIntL.norm hmajL hpt_left
      _ = A * x ^ (σ₁ + 1) * ∫ t in (-T)..(-Tsp), (Real.log (-t + 2)) ^ 2 / t ^ 2 :=
          intervalIntegral.integral_const_mul _ _
      _ = A * x ^ (σ₁ + 1) * ∫ u in Tsp..T, (Real.log (u + 2)) ^ 2 / u ^ 2 := by
          rw [hcomp]
      _ ≤ A * x ^ (σ₁ + 1) * (96 / Real.sqrt Tsp) := by
          gcongr
          exact integral_log_sq_div_sq_le (by linarith) (by linarith)
      _ ≤ A * x ^ (σ₁ + 1) * 96 := by
          have hs1 : (1 : ℝ) ≤ Real.sqrt Tsp := by
            rw [show (1 : ℝ) = Real.sqrt 1 by rw [Real.sqrt_one]]
            exact Real.sqrt_le_sqrt (by linarith)
          gcongr
          rw [div_le_iff₀ (by linarith)]
          nlinarith
      _ = 96 * A * x ^ (σ₁ + 1) := by ring
  have hI2 : ‖∫ t in (-Tsp)..Tsp, F t‖ ≤ (C + L / c₁) * x ^ (σ₁ + 1) * 2 * (2 * Tsp) := by
    have := intervalIntegral.norm_integral_le_of_norm_le_const
      (C := (C + L / c₁) * x ^ (σ₁ + 1) * 2) hpt_mid
    calc ‖∫ t in (-Tsp)..Tsp, F t‖
        ≤ (C + L / c₁) * x ^ (σ₁ + 1) * 2 * |Tsp - (-Tsp)| := this
      _ = (C + L / c₁) * x ^ (σ₁ + 1) * 2 * (2 * Tsp) := by
          rw [show Tsp - (-Tsp) = 2 * Tsp by ring, abs_of_pos (by linarith)]
  -- assemble
  have hLnn : 0 ≤ L := hL.le
  show ‖∫ t in (-T)..T, F t‖ ≤ (192 * A + 4 * Tsp * C + 4 * Tsp / c₁) * (1 + L) * x ^ (σ₁ + 1)
  calc ‖∫ t in (-T)..T, F t‖
      = ‖((∫ t in (-T)..(-Tsp), F t) + ∫ t in (-Tsp)..Tsp, F t) + ∫ t in Tsp..T, F t‖ := by
        rw [hsplit]
    _ ≤ ‖(∫ t in (-T)..(-Tsp), F t) + ∫ t in (-Tsp)..Tsp, F t‖ + ‖∫ t in Tsp..T, F t‖ :=
        norm_add_le _ _
    _ ≤ (‖∫ t in (-T)..(-Tsp), F t‖ + ‖∫ t in (-Tsp)..Tsp, F t‖) + ‖∫ t in Tsp..T, F t‖ := by
        gcongr
        exact norm_add_le _ _
    _ ≤ (96 * A * x ^ (σ₁ + 1) + (C + L / c₁) * x ^ (σ₁ + 1) * 2 * (2 * Tsp))
        + 96 * A * x ^ (σ₁ + 1) := by
        gcongr
    _ = (192 * A + 4 * Tsp * C + 4 * Tsp / c₁ * L) * x ^ (σ₁ + 1) := by ring
    _ ≤ (192 * A + 4 * Tsp * C + 4 * Tsp / c₁) * (1 + L) * x ^ (σ₁ + 1) := by
        apply mul_le_mul_of_nonneg_right _ hxppos.le
        have hTsp0 : (0 : ℝ) ≤ Tsp := by linarith
        have ha : 0 ≤ 192 * A + 4 * Tsp * C :=
          add_nonneg (by positivity) (mul_nonneg (by positivity) hC)
        have hb : 0 ≤ 4 * Tsp / c₁ := by positivity
        have haL : 0 ≤ (192 * A + 4 * Tsp * C) * L := mul_nonneg ha hLnn
        nlinarith [haL, hb]

#print axioms left_line_integral_bound_param
