import TrigPoly5841Full
import SmallTZeroFreeStripT

/-!
# S2947 — clean + UNIFORM 5-8-4-1 zero-free regions (lane 1)

* `riemannZeta_zero_free_5841_clean` : ∃ T₀, ∀ |γ| ≥ T₀, ζ(β+iγ)=0 ⟹ β ≤ 1 − (9/22464)/log|γ|.
  (Absorb 5A₀ ≤ 32·log|γ| via |γ| ≥ exp(5A₀/32): 𝒦 ≤ 864·log|γ|, c₀ = (9/26)/864 = 9/22464 ≈ 4.0e-4 —
  3.6× the forge's 3-4-1 clean constant 1/8960 ≈ 1.1e-4 from the same unoptimized pipeline.)
* `riemannZeta_zero_free_5841_uniform` : ∃ c > 0, ∀ zeros β+iγ, β ≤ 1 − c/log(|γ|+2), combining the clean region with
  the small-|γ| strip (§ SmallTZeroFreeStripT), exactly as `riemannZeta_zero_free_region_uniform` did for 3-4-1.
-/

open Complex

namespace TrigPoly5841

theorem riemannZeta_zero_free_5841_clean :
    ∃ T₀ : ℝ, ∀ γ β : ℝ, T₀ ≤ |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 → β ≤ 1 - (9 / 22464) / Real.log |γ| := by
  obtain ⟨A₀, hA₀nn, hmain⟩ := riemannZeta_zero_free_5841
  refine ⟨max (max (Real.exp 16) (Real.exp (2 * Real.log (16 * (3 + 6 * A₀ + 640)))))
      (Real.exp (5 * A₀ / 32)), ?_⟩
  intro γ β hγ hζ
  have h1 : Real.exp 16 ≤ |γ| := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hγ
  have hexp2 : Real.exp (2 * Real.log (16 * (3 + 6 * A₀ + 640))) ≤ |γ| :=
    le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hγ
  have hexp3 : Real.exp (5 * A₀ / 32) ≤ |γ| := le_trans (le_max_right _ _) hγ
  have h2 : 2 * Real.log (16 * (3 + 6 * A₀ + 640)) ≤ Real.log |γ| := by
    calc 2 * Real.log (16 * (3 + 6 * A₀ + 640))
        = Real.log (Real.exp (2 * Real.log (16 * (3 + 6 * A₀ + 640)))) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp2
  have hL16 : 16 ≤ Real.log |γ| := by
    calc (16 : ℝ) = Real.log (Real.exp 16) := (Real.log_exp 16).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos 16) h1
  have hLpos : 0 < Real.log |γ| := by linarith
  have hL5A : 5 * A₀ / 32 ≤ Real.log |γ| := by
    calc 5 * A₀ / 32 = Real.log (Real.exp (5 * A₀ / 32)) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp3
  have hmain' := hmain γ β h1 h2 hζ
  have hKpos : 0 < 5 * A₀ + 832 * Real.log |γ| := by nlinarith
  have hcmp : 5 * A₀ + 832 * Real.log |γ| ≤ 864 * Real.log |γ| := by nlinarith
  have hdiv : (9 / 26) / (864 * Real.log |γ|) ≤ (9 / 26) / (5 * A₀ + 832 * Real.log |γ|) :=
    div_le_div_of_nonneg_left (by norm_num) hKpos hcmp
  have he : (9 / 22464 : ℝ) / Real.log |γ| = (9 / 26) / (864 * Real.log |γ|) := by
    field_simp; ring
  rw [he]; linarith

theorem riemannZeta_zero_free_5841_uniform :
    ∃ c : ℝ, 0 < c ∧ ∀ β γ : ℝ, riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - c / Real.log (|γ| + 2) := by
  obtain ⟨T₀, hT₀⟩ := riemannZeta_zero_free_5841_clean
  have hc₀ : (0 : ℝ) < 9 / 22464 := by norm_num
  set c₀ : ℝ := 9 / 22464 with hc₀def
  obtain ⟨σ₀, hσhalf, hσ₀, hstrip⟩ := exists_small_t_zero_free_strip_T (max T₀ 2)
    (le_max_of_le_right (by norm_num))
  have hlog2 : (0 : ℝ) < Real.log 2 := Real.log_pos (by norm_num)
  refine ⟨min c₀ ((1 - σ₀) * Real.log 2),
    lt_min hc₀ (mul_pos (by linarith) hlog2), ?_⟩
  intro β γ hζ
  have hre : ((β : ℂ) + I * (γ : ℝ)).re = β := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]
  have him : ((β : ℂ) + I * (γ : ℝ)).im = γ := by
    simp only [Complex.add_im, Complex.ofReal_im, Complex.mul_im, Complex.I_re, Complex.I_im,
      Complex.ofReal_re, one_mul, mul_zero, zero_add]
  have hlogγ2 : Real.log 2 ≤ Real.log (|γ| + 2) :=
    Real.log_le_log (by norm_num) (by linarith [abs_nonneg γ])
  have hlogpos : 0 < Real.log (|γ| + 2) := lt_of_lt_of_le hlog2 hlogγ2
  by_cases hγ : max T₀ 2 ≤ |γ|
  · have h463 := hT₀ γ β (le_trans (le_max_left _ _) hγ) hζ
    have hγ2 : (2 : ℝ) ≤ |γ| := le_trans (le_max_right _ _) hγ
    have hlogγpos : 0 < Real.log |γ| := Real.log_pos (by linarith)
    have hmono : Real.log |γ| ≤ Real.log (|γ| + 2) :=
      Real.log_le_log (by linarith) (by linarith)
    have hfrac : min c₀ ((1 - σ₀) * Real.log 2) / Real.log (|γ| + 2) ≤ c₀ / Real.log |γ| := by
      have h1 : min c₀ ((1 - σ₀) * Real.log 2) / Real.log (|γ| + 2)
          ≤ c₀ / Real.log (|γ| + 2) := by
        gcongr
        exact min_le_left _ _
      have h2 : c₀ / Real.log (|γ| + 2) ≤ c₀ / Real.log |γ| := by
        gcongr
      linarith
    linarith
  · rw [not_le] at hγ
    have hβle : β ≤ σ₀ := by
      by_contra hβ
      rw [not_le] at hβ
      by_cases hβ1 : 1 ≤ β
      · exact riemannZeta_ne_zero_of_one_le_re (by rw [hre]; exact hβ1) hζ
      · rw [not_le] at hβ1
        have hs1 : ((β : ℂ) + I * (γ : ℝ)) ≠ 1 := by
          intro h
          have h' := congrArg Complex.re h
          rw [hre, Complex.one_re] at h'
          linarith
        exact hstrip _ (by rw [hre]; exact hβ) (by rw [hre]; linarith)
          (by rw [him]; linarith) hs1 hζ
    have hc_le : min c₀ ((1 - σ₀) * Real.log 2) / Real.log (|γ| + 2) ≤ 1 - σ₀ := by
      rw [div_le_iff₀ hlogpos]
      calc min c₀ ((1 - σ₀) * Real.log 2) ≤ (1 - σ₀) * Real.log 2 := min_le_right _ _
        _ ≤ (1 - σ₀) * Real.log (|γ| + 2) :=
            mul_le_mul_of_nonneg_left hlogγ2 (by linarith)
    linarith

end TrigPoly5841

#print axioms TrigPoly5841.riemannZeta_zero_free_5841_clean
#print axioms TrigPoly5841.riemannZeta_zero_free_5841_uniform
