import DLVPCharTermNonneg

/-!
# §1467 — DLVP `3-4-1` SERIES positivity for a Dirichlet character (L-analogue of §1358)

Summing the per-term nonnegativity §1466 over `n` (via `Complex.hasSum_re` on the three summable series
`↗Λ`, `↗χ·↗Λ`, `↗χ²·↗Λ`) gives, for `σ > 1`:

  `0 ≤ 3·Re(L(↗Λ) σ) + 4·Re(L(↗χ·↗Λ)(σ+it)) + Re(L(↗χ²·↗Λ)(σ+2it))`,

the `3-4-1` positivity of the (twisted) von Mangoldt L-series = `−ζ'/ζ`, `−L'/L(·,χ)`, `−L'/L(·,χ²)` — the DLVP
`hdlvp` input for the Dirichlet-L zero-free region.  EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction LSeries DirichletCharacter
open scoped LSeries.notation

theorem dlvp_char_series_positivity {N : ℕ} (χ : DirichletCharacter ℂ N) {σ t : ℝ} (hσ : 1 < σ) :
    0 ≤ 3 * (LSeries (↗Λ) ((σ : ℂ))).re
      + 4 * (LSeries (↗χ * ↗Λ) ((σ : ℂ) + I * (t : ℝ))).re
      + (LSeries (↗(χ ^ 2) * ↗Λ) ((σ : ℂ) + I * ((2 * t : ℝ) : ℂ))).re := by
  have hs0 : 1 < ((σ : ℂ)).re := by simpa using hσ
  have hs1 : 1 < ((σ : ℂ) + I * (t : ℝ)).re := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]; simpa using hσ
  have hs2 : 1 < ((σ : ℂ) + I * ((2 * t : ℝ) : ℂ)).re := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]; simpa using hσ
  have H0 := Complex.hasSum_re (LSeriesSummable_vonMangoldt hs0).hasSum
  have H1 := Complex.hasSum_re (LSeriesSummable_twist_vonMangoldt χ hs1).hasSum
  have H2 := Complex.hasSum_re (LSeriesSummable_twist_vonMangoldt (χ ^ 2) hs2).hasSum
  have Hsum := ((H0.mul_left 3).add (H1.mul_left 4)).add H2
  exact Hsum.nonneg (fun n => dlvp_char_term_nonneg χ n)

#print axioms dlvp_char_series_positivity
