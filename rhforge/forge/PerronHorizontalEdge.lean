import PerronFHorizontalBound
import NormLSeriesVonMangoldtMono
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic

set_option maxHeartbeats 1000000

/-!
# §1630 — the horizontal edge integral is `O(1/T²)`: `‖∫_{σ₁}^{σ₀} F(σ'+iT) dσ'‖ ≤ C·(σ₀-σ₁)/T²`

Integrating the per-point horizontal bound §1627 over the fixed-length strip `σ' ∈ [σ₁,σ₀]` (`σ₁>1`), using the
monotonicity §1628 (`‖L(Λ)(σ')‖ ≤ ‖L(Λ)(σ₁)‖`) and `x^{σ'+1} ≤ max(x^{σ₁+1},x^{σ₀+1})`, gives a uniform bound
`‖F(σ'+iT)‖ ≤ ‖L(Λ)(σ₁)‖·max(x^{σ₁+1},x^{σ₀+1})/T²`; `intervalIntegral.norm_integral_le_of_norm_le_const` then
bounds the edge integral by that constant times `|σ₀-σ₁|`.  As `|T|→∞` this `→ 0`: the top/bottom edges of the shift
rectangle vanish, closing the truncated contour.  EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction MeasureTheory
open scoped LSeries.notation

theorem norm_perron_F_horizontal_integral_le {x : ℝ} (hx : 0 < x) {σ₁ σ₀ T : ℝ}
    (hσ₁ : 1 < σ₁) (hσ₀ : σ₁ ≤ σ₀) (hT : T ≠ 0) :
    ‖∫ σ' : ℝ in σ₁..σ₀,
        -logDeriv riemannZeta ((σ' : ℂ) + (T : ℝ) * I) *
            (x : ℂ) ^ (((σ' : ℂ) + (T : ℝ) * I) + 1) /
          (((σ' : ℂ) + (T : ℝ) * I) * (((σ' : ℂ) + (T : ℝ) * I) + 1))‖
      ≤ (‖LSeries (↗Λ) ((σ₁ : ℝ) : ℂ)‖ * max (x ^ (σ₁ + 1)) (x ^ (σ₀ + 1)) / T ^ 2) * |σ₀ - σ₁| := by
  set C : ℝ := ‖LSeries (↗Λ) ((σ₁ : ℝ) : ℂ)‖ * max (x ^ (σ₁ + 1)) (x ^ (σ₀ + 1)) / T ^ 2 with hCdef
  have hT2pos : (0 : ℝ) < T ^ 2 := by positivity
  have hbound : ∀ σ' ∈ Set.Ioc (min σ₁ σ₀) (max σ₁ σ₀),
      ‖-logDeriv riemannZeta ((σ' : ℂ) + (T : ℝ) * I) *
          (x : ℂ) ^ (((σ' : ℂ) + (T : ℝ) * I) + 1) /
        (((σ' : ℂ) + (T : ℝ) * I) * (((σ' : ℂ) + (T : ℝ) * I) + 1))‖ ≤ C := by
    intro σ' hσ'
    rw [min_eq_left hσ₀, max_eq_right hσ₀] at hσ'
    obtain ⟨hlt, hle⟩ := hσ'
    have hσ'1 : 1 < σ' := lt_trans hσ₁ hlt
    -- per-point §1627
    have h27 := norm_perron_F_le_horizontal hx hσ'1 hT
    -- monotonicity of ‖L(Λ)‖
    have hmono : ‖LSeries (↗Λ) ((σ' : ℝ) : ℂ)‖ ≤ ‖LSeries (↗Λ) ((σ₁ : ℝ) : ℂ)‖ :=
      norm_LSeries_vonMangoldt_mono hσ₁ hlt.le
    -- x^{σ'+1} ≤ max
    have hxp : x ^ (σ' + 1) ≤ max (x ^ (σ₁ + 1)) (x ^ (σ₀ + 1)) := by
      rcases le_total 1 x with hx1 | hx1
      · exact le_trans (Real.rpow_le_rpow_of_exponent_le hx1 (by linarith)) (le_max_right _ _)
      · exact le_trans (Real.rpow_le_rpow_of_exponent_ge hx hx1 (by linarith)) (le_max_left _ _)
    have hLnn : (0 : ℝ) ≤ ‖LSeries (↗Λ) ((σ₁ : ℝ) : ℂ)‖ := norm_nonneg _
    have hmaxnn : (0 : ℝ) ≤ max (x ^ (σ₁ + 1)) (x ^ (σ₀ + 1)) :=
      le_trans (Real.rpow_nonneg hx.le _) (le_max_left _ _)
    calc ‖-logDeriv riemannZeta ((σ' : ℂ) + (T : ℝ) * I) *
            (x : ℂ) ^ (((σ' : ℂ) + (T : ℝ) * I) + 1) /
          (((σ' : ℂ) + (T : ℝ) * I) * (((σ' : ℂ) + (T : ℝ) * I) + 1))‖
        ≤ ‖LSeries (↗Λ) ((σ' : ℝ) : ℂ)‖ * x ^ (σ' + 1) / T ^ 2 := h27
      _ ≤ ‖LSeries (↗Λ) ((σ₁ : ℝ) : ℂ)‖ * max (x ^ (σ₁ + 1)) (x ^ (σ₀ + 1)) / T ^ 2 := by
          gcongr
      _ = C := hCdef.symm
  calc ‖∫ σ' : ℝ in σ₁..σ₀,
          -logDeriv riemannZeta ((σ' : ℂ) + (T : ℝ) * I) *
              (x : ℂ) ^ (((σ' : ℂ) + (T : ℝ) * I) + 1) /
            (((σ' : ℂ) + (T : ℝ) * I) * (((σ' : ℂ) + (T : ℝ) * I) + 1))‖
      ≤ C * |σ₀ - σ₁| := intervalIntegral.norm_integral_le_of_norm_le_const hbound

#print axioms norm_perron_F_horizontal_integral_le
