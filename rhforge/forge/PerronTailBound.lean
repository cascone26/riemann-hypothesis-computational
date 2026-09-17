import PerronFVlineIntegrable
import NormLSeriesVonMangoldtLine
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals

set_option maxHeartbeats 4000000

open Complex MeasureTheory Set intervalIntegral

open scoped LSeries.notation ArithmeticFunction

/-!
# §1700 ★★★ — THE PERRON TAIL BOUND: `‖∫_ℝ F − ∫_{|t|≤T} F‖ ≤ 2·‖L(Λ)(σ₀)‖·x^{σ₀+1}/T`

The last analytic estimate of the chain: the two tails `|t| > T` of the full-line Perron integral
are dominated by `‖L(Λ)(σ₀)‖·x^{σ₀+1}/t²` (Dirichlet-series bound on `ζ'/ζ` for `Re = σ₀ > 1`,
`‖s(s+1)‖ ≥ t²`), each integrating to `‖L(Λ)(σ₀)‖·x^{σ₀+1}/T` via `∫_{Ioi T} t^{-2} = 1/T`.
With §1699 (`2πΨ = ∫_ℝ F`) and §1698 (truncated = main + error) this closes the smoothed-PNT
inequality at every finite `T` — optimization over `T` is all that remains.
-/

theorem perron_tail_bound {x σ₀ : ℝ} (hx : 1 ≤ x) (hσ₀ : 1 < σ₀) {T : ℝ} (hT : 1 ≤ T) :
    ‖(∫ t : ℝ,
        -logDeriv riemannZeta ((σ₀ : ℂ) + (t : ℝ) * I) *
            (x : ℂ) ^ (((σ₀ : ℂ) + (t : ℝ) * I) + 1) /
          (((σ₀ : ℂ) + (t : ℝ) * I) * (((σ₀ : ℂ) + (t : ℝ) * I) + 1)))
      - ∫ t in (-T)..T,
        -logDeriv riemannZeta ((σ₀ : ℂ) + (t : ℝ) * I) *
            (x : ℂ) ^ (((σ₀ : ℂ) + (t : ℝ) * I) + 1) /
          (((σ₀ : ℂ) + (t : ℝ) * I) * (((σ₀ : ℂ) + (t : ℝ) * I) + 1))‖
      ≤ 2 * ‖LSeries (↗Λ) ((σ₀ : ℝ) : ℂ)‖ * x ^ (σ₀ + 1) / T := by
  set C : ℝ := ‖LSeries (↗Λ) ((σ₀ : ℝ) : ℂ)‖ with hCdef
  have hC0 : 0 ≤ C := norm_nonneg _
  set F : ℝ → ℂ := fun t =>
    -logDeriv riemannZeta ((σ₀ : ℂ) + (t : ℝ) * I) *
        (x : ℂ) ^ (((σ₀ : ℂ) + (t : ℝ) * I) + 1) /
      (((σ₀ : ℂ) + (t : ℝ) * I) * (((σ₀ : ℂ) + (t : ℝ) * I) + 1)) with hFdef
  have hInt : Integrable F := integrable_perron_F_vline (by linarith) hσ₀
  have hTpos : (0 : ℝ) < T := by linarith
  have hxppos : (0 : ℝ) < x ^ (σ₀ + 1) := Real.rpow_pos_of_pos (by linarith) _
  -- the split: ∫_ℝ − ∫_{-T..T} = ∫_{Iic(-T)} + ∫_{Ioi T}
  have hsplit : (∫ t : ℝ, F t) - (∫ t in (-T)..T, F t)
      = (∫ t in Iic (-T), F t) + ∫ t in Ioi T, F t := by
    have h1 : (∫ t in Iic (-T), F t) + (∫ t in Ioi (-T), F t) = ∫ t : ℝ, F t :=
      integral_Iic_add_Ioi hInt.integrableOn hInt.integrableOn
    have h2 : (∫ t in (-T)..T, F t) + (∫ t in Ioi T, F t) = ∫ t in Ioi (-T), F t :=
      integral_interval_add_Ioi hInt.integrableOn hInt.integrableOn
    linear_combination h1.symm + h2.symm
  rw [hsplit]
  -- the pointwise tail bound for |t| ≥ T
  have hpt : ∀ t : ℝ, T ≤ |t| → ‖F t‖ ≤ C * x ^ (σ₀ + 1) / t ^ 2 := by
    intro t htb
    have htabs : 0 < |t| := lt_of_lt_of_le hTpos htb
    set s : ℂ := (σ₀ : ℂ) + (t : ℝ) * I with hsdef
    have hsre : s.re = σ₀ := by rw [hsdef]; simp
    have hsim : s.im = t := by rw [hsdef]; simp
    have hsre1 : 1 < s.re := by rw [hsre]; exact hσ₀
    -- log-derivative factor
    have hld : ‖logDeriv riemannZeta s‖ ≤ C := by
      rw [logDeriv_riemannZeta_eq_neg_vonMangoldt hsre1, norm_neg]
      have h1 := norm_LSeries_vonMangoldt_le_line hsre1
      rw [hsre] at h1
      exact h1
    -- power factor
    have hxp : ‖(x : ℂ) ^ (s + 1)‖ = x ^ (σ₀ + 1) := by
      rw [Complex.norm_cpow_eq_rpow_re_of_pos (by linarith)]
      congr 1
      rw [Complex.add_re, hsre, Complex.one_re]
    -- denominator
    have hden : t ^ 2 ≤ ‖s * (s + 1)‖ := by
      rw [norm_mul]
      have h1 : |t| ≤ ‖s‖ := by
        calc |t| = |s.im| := by rw [hsim]
          _ ≤ ‖s‖ := Complex.abs_im_le_norm _
      have h2 : |t| ≤ ‖s + 1‖ := by
        have him1 : (s + 1).im = t := by rw [Complex.add_im, hsim, Complex.one_im, add_zero]
        calc |t| = |(s + 1).im| := by rw [him1]
          _ ≤ ‖s + 1‖ := Complex.abs_im_le_norm _
      calc t ^ 2 = |t| * |t| := by rw [← sq_abs]; ring
        _ ≤ ‖s‖ * ‖s + 1‖ := mul_le_mul h1 h2 htabs.le (norm_nonneg _)
    have hdenpos : (0 : ℝ) < ‖s * (s + 1)‖ :=
      lt_of_lt_of_le (pow_two_pos_of_ne_zero (abs_pos.mp htabs)) hden
    rw [hFdef]
    simp only []
    rw [norm_div, norm_mul, norm_neg, hxp]
    have step1 : ‖logDeriv riemannZeta s‖ * x ^ (σ₀ + 1) / ‖s * (s + 1)‖
        ≤ C * x ^ (σ₀ + 1) / ‖s * (s + 1)‖ := by
      gcongr
    have step2 : C * x ^ (σ₀ + 1) / ‖s * (s + 1)‖ ≤ C * x ^ (σ₀ + 1) / t ^ 2 :=
      div_le_div_of_nonneg_left (mul_nonneg hC0 hxppos.le)
        (pow_two_pos_of_ne_zero (abs_pos.mp htabs)) hden
    exact le_trans step1 step2
  -- majorant integrability on Ioi T
  have hmaj : IntegrableOn (fun t : ℝ => C * x ^ (σ₀ + 1) * t ^ (-2 : ℝ)) (Ioi T) := by
    exact (integrableOn_Ioi_rpow_of_lt (by norm_num) hTpos).const_mul _
  have hmajval : (∫ t in Ioi T, C * x ^ (σ₀ + 1) * t ^ (-2 : ℝ))
      = C * x ^ (σ₀ + 1) / T := by
    rw [MeasureTheory.integral_const_mul, integral_Ioi_rpow_of_lt (by norm_num) hTpos]
    rw [show (-2 : ℝ) + 1 = -1 by norm_num, Real.rpow_neg_one]
    field_simp
  have hconv : ∀ t : ℝ, T ≤ t → C * x ^ (σ₀ + 1) / t ^ 2
      = C * x ^ (σ₀ + 1) * t ^ (-2 : ℝ) := by
    intro t ht
    have htpos : 0 < t := lt_of_lt_of_le hTpos ht
    rw [Real.rpow_neg htpos.le, div_eq_mul_inv]
    congr 1
    rw [show (2 : ℝ) = ((2 : ℕ) : ℝ) by norm_num, Real.rpow_natCast]
  -- right tail
  have hIoi : ‖∫ t in Ioi T, F t‖ ≤ C * x ^ (σ₀ + 1) / T := by
    calc ‖∫ t in Ioi T, F t‖ ≤ ∫ t in Ioi T, ‖F t‖ := norm_integral_le_integral_norm _
      _ ≤ ∫ t in Ioi T, C * x ^ (σ₀ + 1) * t ^ (-2 : ℝ) := by
          apply setIntegral_mono_on hInt.norm.integrableOn hmaj measurableSet_Ioi
          intro t ht
          rw [mem_Ioi] at ht
          rw [← hconv t ht.le]
          exact hpt t (by rw [abs_of_pos (lt_of_lt_of_le hTpos ht.le)]; exact ht.le)
      _ = C * x ^ (σ₀ + 1) / T := hmajval
  -- left tail via reflection
  have hIic : ‖∫ t in Iic (-T), F t‖ ≤ C * x ^ (σ₀ + 1) / T := by
    have hrefl : (∫ t in Iic (-T), ‖F t‖) = ∫ u in Ioi T, ‖F (-u)‖ := by
      have := integral_comp_neg_Iic (-T) (fun u : ℝ => ‖F (-u)‖)
      simp only [neg_neg] at this
      rw [← this]
    have hFnegInt : IntegrableOn (fun u : ℝ => ‖F (-u)‖) (Ioi T) := by
      have h1 : Integrable (fun u : ℝ => F (-u)) := by
        have := hInt.comp_neg
        simpa using this
      exact h1.norm.integrableOn
    calc ‖∫ t in Iic (-T), F t‖ ≤ ∫ t in Iic (-T), ‖F t‖ := norm_integral_le_integral_norm _
      _ = ∫ u in Ioi T, ‖F (-u)‖ := hrefl
      _ ≤ ∫ u in Ioi T, C * x ^ (σ₀ + 1) * u ^ (-2 : ℝ) := by
          apply setIntegral_mono_on hFnegInt hmaj measurableSet_Ioi
          intro u hu
          rw [mem_Ioi] at hu
          have hupos : 0 < u := lt_of_lt_of_le hTpos hu.le
          have h1 := hpt (-u) (by rw [abs_neg, abs_of_pos hupos]; exact hu.le)
          rw [show ((-u) : ℝ) ^ 2 = u ^ 2 by ring] at h1
          rw [← hconv u hu.le]
          exact h1
      _ = C * x ^ (σ₀ + 1) / T := hmajval
  calc ‖(∫ t in Iic (-T), F t) + ∫ t in Ioi T, F t‖
      ≤ ‖∫ t in Iic (-T), F t‖ + ‖∫ t in Ioi T, F t‖ := norm_add_le _ _
    _ ≤ C * x ^ (σ₀ + 1) / T + C * x ^ (σ₀ + 1) / T := add_le_add hIic hIoi
    _ = 2 * C * x ^ (σ₀ + 1) / T := by ring

#print axioms perron_tail_bound
