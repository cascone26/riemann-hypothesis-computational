import LogDerivZetaVonMangoldt

/-!
# §1428 — growth bound `Re(-logDeriv ζ(σ+2iγ)) ≤ Re(-logDeriv ζ σ)` (`g₀`-free `hgrowth`)

`-logDeriv ζ = L↗Λ` (§1356) with `Λ ≥ 0`, so termwise `Re(term ↗Λ (σ+2iγ) n) ≤ ‖term ↗Λ (σ+2iγ) n‖ =
‖term ↗Λ σ n‖ = (term ↗Λ σ n).re` (the last since the real-`σ` term is nonnegative real, `LSeries.term_nonneg`).
Summing (`HasSum.le` on the `Complex.hasSum_re` of the von Mangoldt series) gives the growth bound.  Composed with
§1427 (`hpole_relaxed`) this supplies §1426's `hgrowth` with `C := 1/(σ-1)+A`, entirely `g₀`-free.
EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction LSeries
open scoped LSeries.notation ComplexOrder

theorem re_neg_logDeriv_zeta_shift_le {σ γ : ℝ} (hσ : 1 < σ) :
    (-logDeriv riemannZeta ((σ : ℂ) + I * ((2 * γ : ℝ) : ℂ))).re
      ≤ (-logDeriv riemannZeta ((σ : ℂ))).re := by
  have hs0 : 1 < ((σ : ℂ)).re := by simpa using hσ
  have hs2 : 1 < ((σ : ℂ) + I * ((2 * γ : ℝ) : ℂ)).re := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]; simpa using hσ
  have E0 : -logDeriv riemannZeta ((σ : ℂ)) = LSeries (↗Λ) ((σ : ℂ)) := by
    rw [logDeriv_riemannZeta_eq_neg_vonMangoldt hs0]; ring
  have E2 : -logDeriv riemannZeta ((σ : ℂ) + I * ((2 * γ : ℝ) : ℂ))
      = LSeries (↗Λ) ((σ : ℂ) + I * ((2 * γ : ℝ) : ℂ)) := by
    rw [logDeriv_riemannZeta_eq_neg_vonMangoldt hs2]; ring
  rw [E0, E2]
  have H0 := Complex.hasSum_re (LSeriesSummable_vonMangoldt hs0).hasSum
  have H2 := Complex.hasSum_re (LSeriesSummable_vonMangoldt hs2).hasSum
  refine hasSum_le (fun n => ?_) H2 H0
  have hΛnn : (0 : ℂ) ≤ (↗Λ) n := by
    show (0 : ℂ) ≤ ((Λ n : ℝ) : ℂ); exact Complex.zero_le_real.mpr vonMangoldt_nonneg
  have hnn : (0 : ℂ) ≤ LSeries.term (↗Λ) (σ : ℂ) n := LSeries.term_nonneg hΛnn σ
  have hre_eq : (LSeries.term (↗Λ) (σ : ℂ) n).re = ‖LSeries.term (↗Λ) (σ : ℂ) n‖ := by
    conv_lhs => rw [Complex.eq_coe_norm_of_nonneg hnn]
    exact Complex.ofReal_re _
  have hnorm_eq : ‖LSeries.term (↗Λ) ((σ : ℂ) + I * ((2 * γ : ℝ) : ℂ)) n‖
      = ‖LSeries.term (↗Λ) (σ : ℂ) n‖ := by
    rw [LSeries.norm_term_eq, LSeries.norm_term_eq]
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]
  calc (LSeries.term (↗Λ) ((σ : ℂ) + I * ((2 * γ : ℝ) : ℂ)) n).re
      ≤ ‖LSeries.term (↗Λ) ((σ : ℂ) + I * ((2 * γ : ℝ) : ℂ)) n‖ := Complex.re_le_norm _
    _ = ‖LSeries.term (↗Λ) (σ : ℂ) n‖ := hnorm_eq
    _ = (LSeries.term (↗Λ) (σ : ℂ) n).re := hre_eq.symm

#print axioms re_neg_logDeriv_zeta_shift_le
