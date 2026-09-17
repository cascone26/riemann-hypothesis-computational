import NormLSeriesVonMangoldtLine
import LogDerivZetaVonMangoldt
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv
import Mathlib.NumberTheory.LSeries.Deriv
import Mathlib.NumberTheory.LSeries.Dirichlet

set_option maxHeartbeats 1000000

/-!
# §1625 — full-line integrability of the Perron integrand `F` on a vertical line `Re = σ > 1`

`t ↦ F(σ+it) = -ζ'/ζ(σ+it)·x^{(σ+it)+1}/((σ+it)(σ+it+1))` is integrable over all of `ℝ` for `σ > 1`: the numerator
`‖-ζ'/ζ‖ ≤ ‖L(Λ)(σ)‖` (§1591), `‖x^{s+1}‖ = x^{σ+1}` (`norm_cpow_eq_rpow_re_of_pos`), and `‖s‖‖s+1‖ ≥ (1+t²)/2`,
so `‖F(σ+it)‖ ≤ 2‖L(Λ)(σ)‖x^{σ+1}/(1+t²)`, dominated by `integrable_inv_one_add_sq`.  This is the convergence of
the (untruncated) Perron line integral and the source of the `T→∞` vertical-tail decay in the contour shift.
EXIT0 / std-axioms.
-/

open Complex MeasureTheory Filter ArithmeticFunction
open scoped LSeries.notation

theorem integrable_perron_F_vline {x : ℝ} (hx : 0 < x) {σ : ℝ} (hσ : 1 < σ) :
    Integrable (fun t : ℝ =>
      -logDeriv riemannZeta ((σ : ℂ) + (t : ℝ) * I) *
          (x : ℂ) ^ (((σ : ℂ) + (t : ℝ) * I) + 1) /
        (((σ : ℂ) + (t : ℝ) * I) * (((σ : ℂ) + (t : ℝ) * I) + 1))) := by
  set C : ℝ := ‖LSeries (↗Λ) ((σ : ℝ) : ℂ)‖ with hCdef
  have hC0 : 0 ≤ C := norm_nonneg _
  set K : ℝ := 2 * C * x ^ (σ + 1) with hKdef
  have hsre : ∀ t : ℝ, ((σ : ℂ) + (t : ℝ) * I).re = σ := by
    intro t; simp [Complex.add_re, Complex.mul_re]
  have hsne : ∀ t : ℝ, (σ : ℂ) + (t : ℝ) * I ≠ 0 := by
    intro t h; have := congrArg Complex.re h; rw [hsre] at this; simp at this; linarith
  have hs1ne : ∀ t : ℝ, (σ : ℂ) + (t : ℝ) * I + 1 ≠ 0 := by
    intro t h; have := congrArg Complex.re h
    rw [Complex.add_re, hsre, Complex.one_re] at this; simp at this; linarith
  -- continuity ⟹ AEStronglyMeasurable
  have hcont : Continuous (fun t : ℝ =>
      -logDeriv riemannZeta ((σ : ℂ) + (t : ℝ) * I) *
          (x : ℂ) ^ (((σ : ℂ) + (t : ℝ) * I) + 1) /
        (((σ : ℂ) + (t : ℝ) * I) * (((σ : ℂ) + (t : ℝ) * I) + 1))) := by
    have hline : Continuous (fun t : ℝ => (σ : ℂ) + (t : ℝ) * I) := by fun_prop
    have hmaps : Set.MapsTo (fun t : ℝ => (σ : ℂ) + (t : ℝ) * I) Set.univ {s : ℂ | s ≠ 1} := by
      intro t _; show (σ : ℂ) + (t : ℝ) * I ≠ 1
      intro h; have := congrArg Complex.re h; rw [hsre] at this; simp at this; linarith
    have hlogd : Continuous (fun t : ℝ => logDeriv riemannZeta ((σ : ℂ) + (t : ℝ) * I)) := by
      have h1 : Continuous (fun t : ℝ => deriv riemannZeta ((σ : ℂ) + (t : ℝ) * I)) :=
        (analyticOn_riemannZeta.deriv).continuousOn.comp_continuous hline
          (fun t => hmaps (Set.mem_univ t))
      have h2 : Continuous (fun t : ℝ => riemannZeta ((σ : ℂ) + (t : ℝ) * I)) :=
        differentiableOn_riemannZeta.continuousOn.comp_continuous hline
          (fun t => hmaps (Set.mem_univ t))
      have h3 : ∀ t : ℝ, riemannZeta ((σ : ℂ) + (t : ℝ) * I) ≠ 0 := fun t =>
        riemannZeta_ne_zero_of_one_lt_re (by rw [hsre]; exact hσ)
      simp only [logDeriv_apply]; exact h1.div h2 h3
    have hxne : (x : ℂ) ≠ 0 := by exact_mod_cast ne_of_gt hx
    have hcpow : Continuous (fun t : ℝ => (x : ℂ) ^ (((σ : ℂ) + (t : ℝ) * I) + 1)) :=
      (hline.add continuous_const).const_cpow (Or.inl hxne)
    exact (hlogd.neg.mul hcpow).div (hline.mul (hline.add continuous_const))
      (fun t => mul_ne_zero (hsne t) (hs1ne t))
  apply MeasureTheory.Integrable.mono'
    ((integrable_inv_one_add_sq).const_mul K) hcont.aestronglyMeasurable
  refine Filter.Eventually.of_forall (fun t => ?_)
  have h1t : (0 : ℝ) < 1 + t ^ 2 := by positivity
  have hnum : ‖logDeriv riemannZeta ((σ : ℂ) + (t : ℝ) * I)‖ ≤ C := by
    have hsgt : 1 < ((σ : ℂ) + (t : ℝ) * I).re := by rw [hsre]; exact hσ
    rw [logDeriv_riemannZeta_eq_neg_vonMangoldt hsgt, norm_neg]
    have h91 := norm_LSeries_vonMangoldt_le_line hsgt
    rwa [hsre] at h91
  have hcpownorm : ‖(x : ℂ) ^ (((σ : ℂ) + (t : ℝ) * I) + 1)‖ = x ^ (σ + 1) := by
    rw [Complex.norm_cpow_eq_rpow_re_of_pos hx, Complex.add_re, hsre, Complex.one_re]
  have hsy_ge_y : |t| ≤ ‖(σ : ℂ) + (t : ℝ) * I‖ := by
    have := Complex.abs_im_le_norm ((σ : ℂ) + (t : ℝ) * I)
    simpa [Complex.add_im, Complex.mul_im] using this
  have hsy_ge_s : σ ≤ ‖(σ : ℂ) + (t : ℝ) * I‖ := by
    have := Complex.abs_re_le_norm ((σ : ℂ) + (t : ℝ) * I)
    rw [hsre, abs_of_pos (by linarith)] at this; exact this
  have hsy1_ge_y : |t| ≤ ‖(σ : ℂ) + (t : ℝ) * I + 1‖ := by
    have := Complex.abs_im_le_norm ((σ : ℂ) + (t : ℝ) * I + 1)
    simpa [Complex.add_im, Complex.mul_im, Complex.one_im] using this
  have hsy1_ge_s : σ + 1 ≤ ‖(σ : ℂ) + (t : ℝ) * I + 1‖ := by
    have := Complex.abs_re_le_norm ((σ : ℂ) + (t : ℝ) * I + 1)
    rw [Complex.add_re, hsre, Complex.one_re, abs_of_pos (by linarith)] at this; exact this
  have hden_pos : 0 < ‖(σ : ℂ) + (t : ℝ) * I‖ * ‖(σ : ℂ) + (t : ℝ) * I + 1‖ :=
    mul_pos (norm_pos_iff.mpr (hsne t)) (norm_pos_iff.mpr (hs1ne t))
  have hden_ge : (1 + t ^ 2) ≤ 2 * (‖(σ : ℂ) + (t : ℝ) * I‖ * ‖(σ : ℂ) + (t : ℝ) * I + 1‖) := by
    have hyy : |t| * |t| ≤ ‖(σ : ℂ) + (t : ℝ) * I‖ * ‖(σ : ℂ) + (t : ℝ) * I + 1‖ :=
      mul_le_mul hsy_ge_y hsy1_ge_y (abs_nonneg t) (norm_nonneg _)
    have hss : σ * (σ + 1) ≤ ‖(σ : ℂ) + (t : ℝ) * I‖ * ‖(σ : ℂ) + (t : ℝ) * I + 1‖ :=
      mul_le_mul hsy_ge_s hsy1_ge_s (by linarith) (norm_nonneg _)
    have hy2 : |t| * |t| = t ^ 2 := by rw [abs_mul_abs_self, sq]
    rw [hy2] at hyy; nlinarith [hyy, hss]
  have hxpow : (0 : ℝ) ≤ x ^ (σ + 1) := Real.rpow_nonneg hx.le _
  rw [norm_div, norm_mul, norm_mul, norm_neg, hcpownorm, div_le_iff₀ hden_pos]
  have hinvpos : 0 < (1 + t ^ 2)⁻¹ := inv_pos.mpr h1t
  have hstep : (1 : ℝ) ≤ 2 * (1 + t ^ 2)⁻¹ * (‖(σ : ℂ) + (t : ℝ) * I‖ * ‖(σ : ℂ) + (t : ℝ) * I + 1‖) := by
    have h := mul_le_mul_of_nonneg_left hden_ge (le_of_lt hinvpos)
    rw [inv_mul_cancel₀ (ne_of_gt h1t)] at h; nlinarith [h]
  calc ‖logDeriv riemannZeta ((σ : ℂ) + (t : ℝ) * I)‖ * x ^ (σ + 1)
      ≤ C * x ^ (σ + 1) := by
        apply mul_le_mul_of_nonneg_right hnum hxpow
    _ = C * x ^ (σ + 1) * 1 := by ring
    _ ≤ C * x ^ (σ + 1) *
          (2 * (1 + t ^ 2)⁻¹ * (‖(σ : ℂ) + (t : ℝ) * I‖ * ‖(σ : ℂ) + (t : ℝ) * I + 1‖)) := by
        apply mul_le_mul_of_nonneg_left hstep (by positivity)
    _ = K * (1 + t ^ 2)⁻¹ * (‖(σ : ℂ) + (t : ℝ) * I‖ * ‖(σ : ℂ) + (t : ℝ) * I + 1‖) := by
        rw [hKdef]; ring

#print axioms integrable_perron_F_vline
