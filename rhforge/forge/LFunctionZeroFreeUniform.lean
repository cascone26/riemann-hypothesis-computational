import SmallTZeroFreeStripLFunction
import LFunctionZeroFreeAll

set_option maxHeartbeats 1000000

/-!
# §1577 ★★★★★ — THE UNIFORM Dirichlet-L ZERO-FREE REGION: `β ≤ 1 - c/log(|γ|+2)` for EVERY zero of EVERY L

L-analogue of §1575, the family capstone.  For every nontrivial character `χ` mod `N` (quadratic or not),
combining §1568 (`LFunction_zero_free_region_clean_all`, sharp large-`|γ|`) with the strip §1576 (at height
`T = max T₀ 2`) and the 1-line nonvanishing:

  `∃ c > 0, ∀ β γ, LFunction χ(β+iγ) = 0 → β ≤ 1 - c/log(|γ|+2)`   — EVERY zero, NO threshold.

Together with §1575 (ζ), this is the complete classical uniform de la Vallée Poussin zero-free region for the
Riemann zeta function and all Dirichlet L-functions, in the textbook threshold-free shape.  `c` is ineffective
(compactness `σ₀`).  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem LFunction_zero_free_region_uniform {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) :
    ∃ c : ℝ, 0 < c ∧ ∀ β γ : ℝ, LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - c / Real.log (|γ| + 2) := by
  obtain ⟨T₀, c₀, hc₀, hT₀⟩ := LFunction_zero_free_region_clean_all hχ
  obtain ⟨σ₀, hσhalf, hσ₀, hstrip⟩ := exists_small_t_zero_free_strip_LFunction hχ (max T₀ 2)
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
  · have h468 := hT₀ γ β (le_trans (le_max_left _ _) hγ) hζ
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
      · exact LFunction_ne_zero_of_one_le_re χ (Or.inl hχ) (by rw [hre]; exact hβ1) hζ
      · rw [not_le] at hβ1
        exact hstrip _ (by rw [hre]; exact hβ) (by rw [hre]; linarith)
          (by rw [him]; linarith) hζ
    have hc_le : min c₀ ((1 - σ₀) * Real.log 2) / Real.log (|γ| + 2) ≤ 1 - σ₀ := by
      rw [div_le_iff₀ hlogpos]
      calc min c₀ ((1 - σ₀) * Real.log 2) ≤ (1 - σ₀) * Real.log 2 := min_le_right _ _
        _ ≤ (1 - σ₀) * Real.log (|γ| + 2) :=
            mul_le_mul_of_nonneg_left hlogγ2 (by linarith)
    linarith

#print axioms LFunction_zero_free_region_uniform
