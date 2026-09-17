import SmallTZeroFreeStripT
import ZetaZeroFreeClean

set_option maxHeartbeats 1000000

/-!
# §1575 ★★★★★ — THE UNIFORM ζ ZERO-FREE REGION: `β ≤ 1 - c/log(|γ|+2)` for EVERY zero

Combines the sharp large-`|γ|` region §1463 (`β ≤ 1-c₀/log|γ|` for `T₀ ≤ |γ|`) with the small-|t| strip §1574
(at height `T = max T₀ 2`) into the classical UNIFORM statement, valid for **every** zero of ζ with NO threshold:

  `∃ c > 0, ∀ β γ, ζ(β+iγ) = 0 → β ≤ 1 - c/log(|γ|+2)`.

Case split: `|γ| ≥ max T₀ 2` → §1463 (`log` monotone, `c ≤ c₀`); `|γ| < max T₀ 2` → the zero has `β ≤ σ₀` (else
the strip §1574 or the 1-line nonvanishing is contradicted), and `c ≤ (1-σ₀)·log 2` finishes.  Zeros with `β < 1/2`
(incl. the trivial zeros) satisfy the bound trivially since `c/log(|γ|+2) ≤ c/log 2 ≤ 1-σ₀ ≤ 1/2`.
The `σ₀` inherited from §1574 is INEFFECTIVE (compactness), so `c` is ineffective — but the statement is the true
classical uniform de la Vallée Poussin region, threshold-free.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem riemannZeta_zero_free_region_uniform :
    ∃ c : ℝ, 0 < c ∧ ∀ β γ : ℝ, riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - c / Real.log (|γ| + 2) := by
  obtain ⟨T₀, c₀, hc₀, hT₀⟩ := riemannZeta_zero_free_region_clean
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
  · -- large |γ|: sharp region §1463
    have h463 := hT₀ γ β (le_trans (le_max_left _ _) hγ) hζ
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
  · -- small |γ|: β ≤ σ₀
    rw [not_le] at hγ
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

#print axioms riemannZeta_zero_free_region_uniform
