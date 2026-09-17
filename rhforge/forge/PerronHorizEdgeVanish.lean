import PerronHorizontalEdge
import Mathlib.Topology.Order.Bornology

set_option maxHeartbeats 1000000

open Complex ArithmeticFunction Filter Topology
open scoped LSeries.notation

/-- §1637 — the top/bottom rectangle edges VANISH: the Perron horizontal-edge integral → 0 as `T → ∞`.
Squeeze the §1630 bound `≤ (‖L(Λ)(σ₁)‖·max(x^{σ₁+1},x^{σ₀+1})/T²)·|σ₀-σ₁|`, whose right side is `O(1/T²) → 0`.
This is the truncated contour's horizontal edges disappearing — the last ingredient (with the residue and the
vertical lines) of the contour shift. -/
theorem perron_horizontal_edge_tendsto_zero {x σ₁ σ₀ : ℝ} (hx : 0 < x) (hσ₁ : 1 < σ₁) (hσ₀ : σ₁ ≤ σ₀) :
    Tendsto (fun T : ℝ => ∫ σ' : ℝ in σ₁..σ₀,
        -logDeriv riemannZeta ((σ' : ℂ) + (T : ℝ) * I) *
            (x : ℂ) ^ (((σ' : ℂ) + (T : ℝ) * I) + 1) /
          (((σ' : ℂ) + (T : ℝ) * I) * (((σ' : ℂ) + (T : ℝ) * I) + 1)))
      atTop (𝓝 0) := by
  have hg : Tendsto (fun T : ℝ =>
      (‖LSeries (↗Λ) ((σ₁ : ℝ) : ℂ)‖ * max (x ^ (σ₁ + 1)) (x ^ (σ₀ + 1)) / T ^ 2) * |σ₀ - σ₁|)
      atTop (𝓝 0) := by
    have h1 : Tendsto (fun T : ℝ =>
        ‖LSeries (↗Λ) ((σ₁ : ℝ) : ℂ)‖ * max (x ^ (σ₁ + 1)) (x ^ (σ₀ + 1)) / T ^ 2) atTop (𝓝 0) :=
      Tendsto.div_atTop tendsto_const_nhds (tendsto_pow_atTop two_ne_zero)
    simpa using h1.mul_const |σ₀ - σ₁|
  have h : ∀ᶠ T : ℝ in atTop, ‖∫ σ' : ℝ in σ₁..σ₀,
        -logDeriv riemannZeta ((σ' : ℂ) + (T : ℝ) * I) *
            (x : ℂ) ^ (((σ' : ℂ) + (T : ℝ) * I) + 1) /
          (((σ' : ℂ) + (T : ℝ) * I) * (((σ' : ℂ) + (T : ℝ) * I) + 1))‖
      ≤ (‖LSeries (↗Λ) ((σ₁ : ℝ) : ℂ)‖ * max (x ^ (σ₁ + 1)) (x ^ (σ₀ + 1)) / T ^ 2) * |σ₀ - σ₁| := by
    filter_upwards [eventually_ne_atTop 0] with T hT
    exact norm_perron_F_horizontal_integral_le hx hσ₁ hσ₀ hT
  exact squeeze_zero_norm' h hg

#print axioms perron_horizontal_edge_tendsto_zero
