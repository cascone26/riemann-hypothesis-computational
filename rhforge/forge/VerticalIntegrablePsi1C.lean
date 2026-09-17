import Psi1CMellinForm
import NormLSeriesVonMangoldtLine
import LogDerivZetaVonMangoldt
import Mathlib.Analysis.MellinTransform
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals
import Mathlib.NumberTheory.LSeries.Deriv

set_option maxHeartbeats 1000000

/-!
# §1600 — `VerticalIntegrable (mellin psi1C) (-σ₀-1)` for `σ₀ > 1`  (the meaty inversion hypothesis)

On the line `Re = -σ₀-1`, §1598 gives `mellin Ψ = -ζ'/ζ(s_y)/(s_y(s_y+1))` with `s_y = σ₀ - iy`.  The numerator is
bounded by `‖L(Λ)(σ₀)‖` (§1591 line bound + §1356), and `‖s_y‖·‖s_y+1‖ ≥ (1+y²)/2` (each factor `≥ |y|` and `≥ σ₀`,
`σ₀(σ₀+1) > 2`).  So `‖mellin Ψ‖ ≤ 2‖L(Λ)(σ₀)‖/(1+y²)`, integrable (`integrable_inv_one_add_sq`).  Measurability
comes from the explicit formula (continuous: `analyticOn_riemannZeta.deriv` + nonvanishing).  EXIT0 / std-axioms.
-/

open Complex Chebyshev MeasureTheory Filter ArithmeticFunction
open scoped LSeries.notation

theorem verticalIntegrable_mellin_psi1C {σ₀ : ℝ} (hσ₀ : 1 < σ₀) :
    Complex.VerticalIntegrable (mellin psi1C) (-σ₀ - 1) := by
  rw [Complex.VerticalIntegrable]
  set C : ℝ := ‖LSeries (↗Λ) ((σ₀ : ℝ) : ℂ)‖ with hCdef
  have hC0 : 0 ≤ C := norm_nonneg _
  -- the explicit formula G y
  set G : ℝ → ℂ := fun y => -logDeriv riemannZeta ((σ₀ : ℂ) - (y : ℝ) * I)
      / (((σ₀ : ℂ) - (y : ℝ) * I) * (((σ₀ : ℂ) - (y : ℝ) * I) + 1)) with hGdef
  have hsyre : ∀ y : ℝ, ((σ₀ : ℂ) - (y : ℝ) * I).re = σ₀ := by
    intro y
    simp [Complex.sub_re, Complex.mul_re]
  have hsyne : ∀ y : ℝ, (σ₀ : ℂ) - (y : ℝ) * I ≠ 0 := by
    intro y h
    have := congrArg Complex.re h
    rw [hsyre] at this; simp at this; linarith
  have hsy1ne : ∀ y : ℝ, (σ₀ : ℂ) - (y : ℝ) * I + 1 ≠ 0 := by
    intro y h
    have := congrArg Complex.re h
    rw [Complex.add_re, hsyre, Complex.one_re] at this; simp at this; linarith
  -- integrand equals G
  have hEq : (fun y : ℝ => mellin psi1C (((-σ₀ - 1 : ℝ) : ℂ) + (y : ℝ) * I)) = G := by
    funext y
    have hsygt : 1 < ((σ₀ : ℂ) - (y : ℝ) * I).re := by rw [hsyre]; exact hσ₀
    have harg : ((-σ₀ - 1 : ℝ) : ℂ) + (y : ℝ) * I = -((σ₀ : ℂ) - (y : ℝ) * I) - 1 := by
      push_cast; ring
    rw [harg]
    have h98 := mellin_psi1C_eq hsygt
    rw [hGdef]
    rw [eq_div_iff (mul_ne_zero (hsyne y) (hsy1ne y))]
    linear_combination h98
  rw [hEq]
  -- continuity of G  ⟹ AEStronglyMeasurable
  have hcont : Continuous G := by
    have hline : Continuous (fun y : ℝ => (σ₀ : ℂ) - (y : ℝ) * I) := by fun_prop
    have hmaps : Set.MapsTo (fun y : ℝ => (σ₀ : ℂ) - (y : ℝ) * I) Set.univ {s : ℂ | s ≠ 1} := by
      intro y _
      show (σ₀ : ℂ) - (y : ℝ) * I ≠ 1
      intro h; have := congrArg Complex.re h; rw [hsyre] at this
      simp at this; linarith
    have hderiv : ContinuousOn (deriv riemannZeta) {s : ℂ | s ≠ 1} :=
      (analyticOn_riemannZeta.deriv).continuousOn
    have hzeta : ContinuousOn riemannZeta {s : ℂ | s ≠ 1} :=
      differentiableOn_riemannZeta.continuousOn
    have hlogd : Continuous (fun y : ℝ => logDeriv riemannZeta ((σ₀ : ℂ) - (y : ℝ) * I)) := by
      have h1 : Continuous (fun y : ℝ => deriv riemannZeta ((σ₀ : ℂ) - (y : ℝ) * I)) :=
        hderiv.comp_continuous hline (fun y => hmaps (Set.mem_univ y))
      have h2 : Continuous (fun y : ℝ => riemannZeta ((σ₀ : ℂ) - (y : ℝ) * I)) :=
        hzeta.comp_continuous hline (fun y => hmaps (Set.mem_univ y))
      have h3 : ∀ y : ℝ, riemannZeta ((σ₀ : ℂ) - (y : ℝ) * I) ≠ 0 := fun y =>
        riemannZeta_ne_zero_of_one_lt_re (by rw [hsyre]; exact hσ₀)
      simp only [logDeriv_apply]
      exact h1.div h2 h3
    rw [hGdef]
    exact (hlogd.neg).div (hline.mul (hline.add continuous_const))
      (fun y => mul_ne_zero (hsyne y) (hsy1ne y))
  -- domination
  apply MeasureTheory.Integrable.mono'
    ((integrable_inv_one_add_sq).const_mul (2 * C)) hcont.aestronglyMeasurable
  refine Filter.Eventually.of_forall (fun y => ?_)
  have h1y : (0 : ℝ) < 1 + y ^ 2 := by positivity
  have hden_pos : 0 < ‖(σ₀ : ℂ) - (y : ℝ) * I‖ * ‖(σ₀ : ℂ) - (y : ℝ) * I + 1‖ :=
    mul_pos (norm_pos_iff.mpr (hsyne y)) (norm_pos_iff.mpr (hsy1ne y))
  have hnum : ‖logDeriv riemannZeta ((σ₀ : ℂ) - (y : ℝ) * I)‖ ≤ C := by
    have hsygt : 1 < ((σ₀ : ℂ) - (y : ℝ) * I).re := by rw [hsyre]; exact hσ₀
    rw [logDeriv_riemannZeta_eq_neg_vonMangoldt hsygt, norm_neg]
    have h91 := norm_LSeries_vonMangoldt_le_line hsygt
    rwa [hsyre] at h91
  have hsy_ge_y : |y| ≤ ‖(σ₀ : ℂ) - (y : ℝ) * I‖ := by
    have := Complex.abs_im_le_norm ((σ₀ : ℂ) - (y : ℝ) * I)
    simpa [Complex.sub_im, Complex.mul_im] using this
  have hsy_ge_s : σ₀ ≤ ‖(σ₀ : ℂ) - (y : ℝ) * I‖ := by
    have := Complex.abs_re_le_norm ((σ₀ : ℂ) - (y : ℝ) * I)
    rw [hsyre, abs_of_pos (by linarith)] at this; exact this
  have hsy1_ge_y : |y| ≤ ‖(σ₀ : ℂ) - (y : ℝ) * I + 1‖ := by
    have := Complex.abs_im_le_norm ((σ₀ : ℂ) - (y : ℝ) * I + 1)
    simpa [Complex.add_im, Complex.sub_im, Complex.mul_im, Complex.one_im] using this
  have hsy1_ge_s : σ₀ + 1 ≤ ‖(σ₀ : ℂ) - (y : ℝ) * I + 1‖ := by
    have := Complex.abs_re_le_norm ((σ₀ : ℂ) - (y : ℝ) * I + 1)
    rw [Complex.add_re, hsyre, Complex.one_re, abs_of_pos (by linarith)] at this; exact this
  have hden_ge : (1 + y ^ 2) ≤ 2 * (‖(σ₀ : ℂ) - (y : ℝ) * I‖ * ‖(σ₀ : ℂ) - (y : ℝ) * I + 1‖) := by
    have hyy : |y| * |y| ≤ ‖(σ₀ : ℂ) - (y : ℝ) * I‖ * ‖(σ₀ : ℂ) - (y : ℝ) * I + 1‖ :=
      mul_le_mul hsy_ge_y hsy1_ge_y (abs_nonneg y) (norm_nonneg _)
    have hss : σ₀ * (σ₀ + 1) ≤ ‖(σ₀ : ℂ) - (y : ℝ) * I‖ * ‖(σ₀ : ℂ) - (y : ℝ) * I + 1‖ :=
      mul_le_mul hsy_ge_s hsy1_ge_s (by linarith) (norm_nonneg _)
    have hy2 : |y| * |y| = y ^ 2 := by rw [abs_mul_abs_self, sq]
    rw [hy2] at hyy
    nlinarith [hyy, hss]
  have hinvpos : 0 < (1 + y ^ 2)⁻¹ := inv_pos.mpr h1y
  have hinv : (1 + y ^ 2)⁻¹ * (1 + y ^ 2) = 1 := inv_mul_cancel₀ (ne_of_gt h1y)
  have hstep1 : (1 : ℝ) ≤ 2 * (1 + y ^ 2)⁻¹ * (‖(σ₀ : ℂ) - (y : ℝ) * I‖ * ‖(σ₀ : ℂ) - (y : ℝ) * I + 1‖) := by
    have h := mul_le_mul_of_nonneg_left hden_ge (le_of_lt hinvpos)
    rw [hinv] at h
    nlinarith [h]
  rw [hGdef, norm_div, norm_neg, norm_mul, div_le_iff₀ hden_pos]
  calc ‖logDeriv riemannZeta ((σ₀ : ℂ) - (y : ℝ) * I)‖
      ≤ C * (2 * (1 + y ^ 2)⁻¹ * (‖(σ₀ : ℂ) - (y : ℝ) * I‖ * ‖(σ₀ : ℂ) - (y : ℝ) * I + 1‖)) :=
        le_trans hnum (le_mul_of_one_le_right hC0 hstep1)
    _ = 2 * C * (1 + y ^ 2)⁻¹ * (‖(σ₀ : ℂ) - (y : ℝ) * I‖ * ‖(σ₀ : ℂ) - (y : ℝ) * I + 1‖) := by ring

#print axioms verticalIntegrable_mellin_psi1C
