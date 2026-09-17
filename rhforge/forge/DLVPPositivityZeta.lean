import DLVPTermNonneg
import LogDerivZetaVonMangoldt

/-!
# §1358 — DLVP `3-4-1` positivity for `-logDeriv ζ` (§1355's `hdlvp`)

Summing the termwise nonnegativity §1357 over `n` (via `Complex.hasSum_re` + `LSeriesSummable_vonMangoldt`,
combined with `HasSum.mul_left`/`HasSum.add`/`HasSum.nonneg`) and rewriting `-logDeriv ζ = L ↗Λ` (§1356) yields the
classical de la Vallée Poussin positivity for `σ > 1`:

  `0 ≤ 3·Re(-logDeriv ζ σ) + 4·Re(-logDeriv ζ (σ+iγ)) + Re(-logDeriv ζ (σ+2iγ))`.

This is exactly §1355's `hdlvp` (with `T₀ = Re(-logDeriv ζ σ)`, `T₂ = Re(-logDeriv ζ (σ+2iγ))`).  EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction LSeries
open scoped LSeries.notation

theorem dlvp_positivity_riemannZeta {σ γ : ℝ} (hσ : 1 < σ) :
    0 ≤ 3 * (-logDeriv riemannZeta ((σ : ℂ))).re
      + 4 * (-logDeriv riemannZeta ((σ : ℂ) + I * (γ : ℝ))).re
      + (-logDeriv riemannZeta ((σ : ℂ) + I * ((2 * γ : ℝ) : ℂ))).re := by
  have hs0 : 1 < ((σ : ℂ)).re := by simpa using hσ
  have hs1 : 1 < ((σ : ℂ) + I * (γ : ℝ)).re := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]; simpa using hσ
  have hs2 : 1 < ((σ : ℂ) + I * ((2 * γ : ℝ) : ℂ)).re := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]; simpa using hσ
  have H0 := Complex.hasSum_re (LSeriesSummable_vonMangoldt hs0).hasSum
  have H1 := Complex.hasSum_re (LSeriesSummable_vonMangoldt hs1).hasSum
  have H2 := Complex.hasSum_re (LSeriesSummable_vonMangoldt hs2).hasSum
  have Hsum := ((H0.mul_left 3).add (H1.mul_left 4)).add H2
  have hnn := Hsum.nonneg (fun n => dlvp_term_combination_nonneg n)
  have E0 : -logDeriv riemannZeta ((σ : ℂ)) = LSeries (↗Λ) ((σ : ℂ)) := by
    rw [logDeriv_riemannZeta_eq_neg_vonMangoldt hs0]; ring
  have E1 : -logDeriv riemannZeta ((σ : ℂ) + I * (γ : ℝ)) = LSeries (↗Λ) ((σ : ℂ) + I * (γ : ℝ)) := by
    rw [logDeriv_riemannZeta_eq_neg_vonMangoldt hs1]; ring
  have E2 : -logDeriv riemannZeta ((σ : ℂ) + I * ((2 * γ : ℝ) : ℂ))
      = LSeries (↗Λ) ((σ : ℂ) + I * ((2 * γ : ℝ) : ℂ)) := by
    rw [logDeriv_riemannZeta_eq_neg_vonMangoldt hs2]; ring
  rw [E0, E1, E2]
  exact hnn

#print axioms dlvp_positivity_riemannZeta
