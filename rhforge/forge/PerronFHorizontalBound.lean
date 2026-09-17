import NormLSeriesVonMangoldtLine
import LogDerivZetaVonMangoldt
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv
import Mathlib.NumberTheory.LSeries.Deriv

set_option maxHeartbeats 1000000

/-!
# §1627 — horizontal decay bound: `‖F(σ+iT)‖ ≤ ‖L(Λ)(σ)‖·x^{σ+1}/T²`

At height `T` (`T ≠ 0`), the Perron integrand satisfies `‖F(σ+iT)‖ ≤ ‖L(Λ)(σ)‖·x^{σ+1}/T²`: the numerator is
`‖-ζ'/ζ‖·‖x^{s+1}‖ ≤ ‖L(Λ)(σ)‖·x^{σ+1}` (§1591 + `norm_cpow`) and the denominator
`‖(σ+iT)(σ+iT+1)‖ ≥ |T|·|T| = T²` (both factors have imaginary part `T`).  As `|T|→∞` this forces the top/bottom
horizontal edges of the shift rectangle to vanish — the input that closes the truncated contour.  EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction
open scoped LSeries.notation

theorem norm_perron_F_le_horizontal {x : ℝ} (hx : 0 < x) {σ T : ℝ} (hσ : 1 < σ) (hT : T ≠ 0) :
    ‖-logDeriv riemannZeta ((σ : ℂ) + (T : ℝ) * I) *
        (x : ℂ) ^ (((σ : ℂ) + (T : ℝ) * I) + 1) /
      (((σ : ℂ) + (T : ℝ) * I) * (((σ : ℂ) + (T : ℝ) * I) + 1))‖
      ≤ ‖LSeries (↗Λ) ((σ : ℝ) : ℂ)‖ * x ^ (σ + 1) / T ^ 2 := by
  have hsre : ((σ : ℂ) + (T : ℝ) * I).re = σ := by simp [Complex.add_re, Complex.mul_re]
  have hsgt : 1 < ((σ : ℂ) + (T : ℝ) * I).re := by rw [hsre]; exact hσ
  set C : ℝ := ‖LSeries (↗Λ) ((σ : ℝ) : ℂ)‖ with hCdef
  have hC0 : 0 ≤ C := norm_nonneg _
  have hxpow : (0 : ℝ) ≤ x ^ (σ + 1) := Real.rpow_nonneg hx.le _
  -- numerator bound
  have hnum : ‖logDeriv riemannZeta ((σ : ℂ) + (T : ℝ) * I)‖ ≤ C := by
    rw [logDeriv_riemannZeta_eq_neg_vonMangoldt hsgt, norm_neg]
    have h91 := norm_LSeries_vonMangoldt_le_line hsgt
    rwa [hsre] at h91
  have hcpownorm : ‖(x : ℂ) ^ (((σ : ℂ) + (T : ℝ) * I) + 1)‖ = x ^ (σ + 1) := by
    rw [Complex.norm_cpow_eq_rpow_re_of_pos hx, Complex.add_re, hsre, Complex.one_re]
  -- denominator bounds
  have hsne : (σ : ℂ) + (T : ℝ) * I ≠ 0 := by
    intro h; have := congrArg Complex.re h; rw [hsre] at this; simp at this; linarith
  have hs1ne : (σ : ℂ) + (T : ℝ) * I + 1 ≠ 0 := by
    intro h; have := congrArg Complex.re h
    rw [Complex.add_re, hsre, Complex.one_re] at this; simp at this; linarith
  have hs_ge : |T| ≤ ‖(σ : ℂ) + (T : ℝ) * I‖ := by
    have := Complex.abs_im_le_norm ((σ : ℂ) + (T : ℝ) * I)
    simpa [Complex.add_im, Complex.mul_im] using this
  have hs1_ge : |T| ≤ ‖(σ : ℂ) + (T : ℝ) * I + 1‖ := by
    have := Complex.abs_im_le_norm ((σ : ℂ) + (T : ℝ) * I + 1)
    simpa [Complex.add_im, Complex.mul_im, Complex.one_im] using this
  have hTpos : (0 : ℝ) < |T| := abs_pos.mpr hT
  have hden_pos : 0 < ‖(σ : ℂ) + (T : ℝ) * I‖ * ‖(σ : ℂ) + (T : ℝ) * I + 1‖ :=
    mul_pos (norm_pos_iff.mpr hsne) (norm_pos_iff.mpr hs1ne)
  have hT2 : T ^ 2 ≤ ‖(σ : ℂ) + (T : ℝ) * I‖ * ‖(σ : ℂ) + (T : ℝ) * I + 1‖ := by
    have h := mul_le_mul hs_ge hs1_ge (abs_nonneg T) (norm_nonneg _)
    have : |T| * |T| = T ^ 2 := by rw [abs_mul_abs_self, sq]
    rw [this] at h; exact h
  have hT2pos : (0 : ℝ) < T ^ 2 := by positivity
  rw [norm_div, norm_mul, norm_mul, norm_neg, hcpownorm, div_le_div_iff₀ hden_pos hT2pos]
  calc ‖logDeriv riemannZeta ((σ : ℂ) + (T : ℝ) * I)‖ * x ^ (σ + 1) * T ^ 2
      ≤ C * x ^ (σ + 1) * (‖(σ : ℂ) + (T : ℝ) * I‖ * ‖(σ : ℂ) + (T : ℝ) * I + 1‖) := by
        apply mul_le_mul _ hT2 (le_of_lt hT2pos) (by positivity)
        exact mul_le_mul_of_nonneg_right hnum hxpow
    _ = C * x ^ (σ + 1) * (‖(σ : ℂ) + (T : ℝ) * I‖ * ‖(σ : ℂ) + (T : ℝ) * I + 1‖) := rfl

#print axioms norm_perron_F_le_horizontal
