import ZetaZeroFreeClean

/-!
# §1464 — the zero-free region in NAMED nonvanishing form

The standard way the de la Vallée Poussin region is stated and used: `ζ(s) ≠ 0` throughout the region
`Re s > 1 - c₀/log|Im s|` (for `|Im s|` beyond a threshold).  Contrapositive of §1463.  This is the form that
feeds PNT-with-error and every downstream application.  EXIT0 / std-axioms.
-/

open Complex

theorem riemannZeta_ne_zero_in_region :
    ∃ T₀ c₀ : ℝ, 0 < c₀ ∧ ∀ s : ℂ, T₀ ≤ |s.im| → 1 - c₀ / Real.log |s.im| < s.re →
      riemannZeta s ≠ 0 := by
  obtain ⟨T₀, c₀, hc₀, hmain⟩ := riemannZeta_zero_free_region_clean
  refine ⟨T₀, c₀, hc₀, ?_⟩
  intro s hT hσ hζ
  have hs_eq : s = ((s.re : ℝ) : ℂ) + I * ((s.im : ℝ) : ℂ) := by
    rw [mul_comm]; exact (Complex.re_add_im s).symm
  rw [hs_eq] at hζ
  have hle := hmain s.im s.re hT hζ
  linarith

#print axioms riemannZeta_ne_zero_in_region
