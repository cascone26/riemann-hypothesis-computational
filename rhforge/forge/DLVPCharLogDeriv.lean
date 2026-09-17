import DLVPCharSeriesPos

/-!
# §1468 — DLVP `3-4-1` positivity in `logDeriv` form (`−ζ'/ζ`, `−L'/L(·,χ)`, `−L'/L(·,χ²)`)

Rewriting the twisted von Mangoldt L-series of §1467 as logarithmic derivatives via Mathlib's
`LSeries_vonMangoldt_eq_deriv_riemannZeta_div` (`L(↗Λ) s = −ζ'(s)/ζ(s)`) and
`LSeries_twist_vonMangoldt_eq` (`L(↗χ·↗Λ) s = −L'(s,χ)/L(s,χ)`) turns §1467 into the standard analytic form

  `0 ≤ 3·Re(−logDeriv ζ σ) + 4·Re(−logDeriv L(·,χ) (σ+it)) + Re(−logDeriv L(·,χ²) (σ+2it))`,

i.e. the de la Vallée Poussin `hdlvp` input, now directly in terms of `−ζ'/ζ` and `−L'/L`.  This is the exact
brick the Borel–Carathéodory / region tower consumes.  EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction LSeries DirichletCharacter
open scoped LSeries.notation

theorem dlvp_char_positivity_logDeriv {N : ℕ} (χ : DirichletCharacter ℂ N) {σ t : ℝ} (hσ : 1 < σ) :
    0 ≤ 3 * (-logDeriv riemannZeta ((σ : ℂ))).re
      + 4 * (-logDeriv (LSeries (↗χ)) ((σ : ℂ) + I * (t : ℝ))).re
      + (-logDeriv (LSeries (↗(χ ^ 2))) ((σ : ℂ) + I * ((2 * t : ℝ) : ℂ))).re := by
  have hs0 : 1 < ((σ : ℂ)).re := by simpa using hσ
  have hs1 : 1 < ((σ : ℂ) + I * (t : ℝ)).re := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]; simpa using hσ
  have hs2 : 1 < ((σ : ℂ) + I * ((2 * t : ℝ) : ℂ)).re := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]; simpa using hσ
  have E0 : LSeries (↗Λ) ((σ : ℂ)) = -logDeriv riemannZeta ((σ : ℂ)) := by
    rw [LSeries_vonMangoldt_eq_deriv_riemannZeta_div hs0, logDeriv_apply, neg_div]
  have E1 : LSeries (↗χ * ↗Λ) ((σ : ℂ) + I * (t : ℝ))
      = -logDeriv (LSeries (↗χ)) ((σ : ℂ) + I * (t : ℝ)) := by
    rw [LSeries_twist_vonMangoldt_eq χ hs1, logDeriv_apply, neg_div]
  have E2 : LSeries (↗(χ ^ 2) * ↗Λ) ((σ : ℂ) + I * ((2 * t : ℝ) : ℂ))
      = -logDeriv (LSeries (↗(χ ^ 2))) ((σ : ℂ) + I * ((2 * t : ℝ) : ℂ)) := by
    rw [LSeries_twist_vonMangoldt_eq (χ ^ 2) hs2, logDeriv_apply, neg_div]
  rw [← E0, ← E1, ← E2]
  exact dlvp_char_series_positivity χ hσ

#print axioms dlvp_char_positivity_logDeriv
