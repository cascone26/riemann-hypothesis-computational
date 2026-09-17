import NormSumIocTwistUniform

open Finset

/-! §1524 uniform twisted block bound on {δ ≤ Re s, ‖s‖ ≤ R}. EXIT0 / std-axioms. -/

theorem norm_sum_Ioc_twist_unif_region {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M} (hχ : χ ≠ 1)
    {s : ℂ} {δ R : ℝ} (hδ : 0 < δ) (hσ : δ ≤ s.re) (hsR : ‖s‖ ≤ R) {m n : ℕ} (hm : 1 ≤ m) (hmn : m < n) :
    ‖∑ i ∈ Finset.Ioc m n, (i : ℂ) ^ (-s) • χ (i : ZMod M)‖
      ≤ ((n : ℝ) ^ (-δ) + ((m + 1 : ℕ) : ℝ) ^ (-δ) + R * (m : ℝ) ^ (-δ) / δ) * M := by
  have hσpos : 0 < s.re := lt_of_lt_of_le hδ hσ
  have hs0 : s ≠ 0 := by rintro rfl; rw [Complex.zero_re] at hσpos; exact lt_irrefl 0 hσpos
  have h := norm_sum_Ioc_twist_uniform_le hχ hσpos hs0 hm hmn
  have hMnn : (0 : ℝ) ≤ (M : ℝ) := Nat.cast_nonneg M
  have hn1 : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast (by omega : 1 ≤ n)
  have hm1 : (1 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hm
  have hm11 : (1 : ℝ) ≤ ((m + 1 : ℕ) : ℝ) := by exact_mod_cast (by omega : 1 ≤ m + 1)
  have hsRnn : (0 : ℝ) ≤ R := le_trans (norm_nonneg s) hsR
  have e1 : (n : ℝ) ^ (-s.re) ≤ (n : ℝ) ^ (-δ) :=
    Real.rpow_le_rpow_of_exponent_le hn1 (by linarith)
  have e2 : ((m + 1 : ℕ) : ℝ) ^ (-s.re) ≤ ((m + 1 : ℕ) : ℝ) ^ (-δ) :=
    Real.rpow_le_rpow_of_exponent_le hm11 (by linarith)
  have e3 : (m : ℝ) ^ (-s.re) ≤ (m : ℝ) ^ (-δ) :=
    Real.rpow_le_rpow_of_exponent_le hm1 (by linarith)
  have hA : ‖s‖ * (m : ℝ) ^ (-s.re) ≤ R * (m : ℝ) ^ (-δ) :=
    mul_le_mul hsR e3 (Real.rpow_nonneg (by linarith) _) hsRnn
  have htail : ‖s‖ * (m : ℝ) ^ (-s.re) / s.re ≤ R * (m : ℝ) ^ (-δ) / δ := by
    rw [div_eq_mul_inv, div_eq_mul_inv]
    exact mul_le_mul hA (inv_anti₀ hδ hσ) (inv_nonneg.mpr hσpos.le)
      (mul_nonneg hsRnn (Real.rpow_nonneg (by linarith) _))
  calc ‖∑ i ∈ Finset.Ioc m n, (i : ℂ) ^ (-s) • χ (i : ZMod M)‖
      ≤ (n : ℝ) ^ (-s.re) * M + ((m + 1 : ℕ) : ℝ) ^ (-s.re) * M
        + ‖s‖ * (m : ℝ) ^ (-s.re) / s.re * M := h
    _ = ((n : ℝ) ^ (-s.re) + ((m + 1 : ℕ) : ℝ) ^ (-s.re) + ‖s‖ * (m : ℝ) ^ (-s.re) / s.re) * M := by
        ring
    _ ≤ ((n : ℝ) ^ (-δ) + ((m + 1 : ℕ) : ℝ) ^ (-δ) + R * (m : ℝ) ^ (-δ) / δ) * M := by
        apply mul_le_mul_of_nonneg_right _ hMnn
        linarith [e1, e2, htail]

#print axioms norm_sum_Ioc_twist_unif_region
