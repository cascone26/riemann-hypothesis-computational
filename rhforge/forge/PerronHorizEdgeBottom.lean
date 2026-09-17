import PerronHorizontalEdge
import Mathlib.Topology.Order.Bornology

set_option maxHeartbeats 1000000

open Complex ArithmeticFunction Filter Topology
open scoped LSeries.notation

/-- §1646 — the BOTTOM horizontal edge (at height `-T`) also vanishes as `T → ∞`:
`∫_{σ₁}^{σ₀} F(σ'+i(-T)) dσ' → 0`.  Same §1630 squeeze as §1637 but at `Im = -T` (bound uses `(-T)² = T²`).
The bottom side of the truncated shift rectangle disappears — completing (with §1637 top) both horizontal edges. -/
theorem perron_horizontal_edge_bottom_tendsto_zero {x σ₁ σ₀ : ℝ} (hx : 0 < x) (hσ₁ : 1 < σ₁)
    (hσ₀ : σ₁ ≤ σ₀) :
    Tendsto (fun T : ℝ => ∫ σ' : ℝ in σ₁..σ₀,
        -logDeriv riemannZeta ((σ' : ℂ) + (-T : ℝ) * I) *
            (x : ℂ) ^ (((σ' : ℂ) + (-T : ℝ) * I) + 1) /
          (((σ' : ℂ) + (-T : ℝ) * I) * (((σ' : ℂ) + (-T : ℝ) * I) + 1)))
      atTop (𝓝 0) := by
  have hg : Tendsto (fun T : ℝ =>
      (‖LSeries (↗Λ) ((σ₁ : ℝ) : ℂ)‖ * max (x ^ (σ₁ + 1)) (x ^ (σ₀ + 1)) / T ^ 2) * |σ₀ - σ₁|)
      atTop (𝓝 0) := by
    have h1 : Tendsto (fun T : ℝ =>
        ‖LSeries (↗Λ) ((σ₁ : ℝ) : ℂ)‖ * max (x ^ (σ₁ + 1)) (x ^ (σ₀ + 1)) / T ^ 2) atTop (𝓝 0) :=
      Tendsto.div_atTop tendsto_const_nhds (tendsto_pow_atTop two_ne_zero)
    simpa using h1.mul_const |σ₀ - σ₁|
  have h : ∀ᶠ T : ℝ in atTop, ‖∫ σ' : ℝ in σ₁..σ₀,
        -logDeriv riemannZeta ((σ' : ℂ) + (-T : ℝ) * I) *
            (x : ℂ) ^ (((σ' : ℂ) + (-T : ℝ) * I) + 1) /
          (((σ' : ℂ) + (-T : ℝ) * I) * (((σ' : ℂ) + (-T : ℝ) * I) + 1))‖
      ≤ (‖LSeries (↗Λ) ((σ₁ : ℝ) : ℂ)‖ * max (x ^ (σ₁ + 1)) (x ^ (σ₀ + 1)) / T ^ 2) * |σ₀ - σ₁| := by
    filter_upwards [eventually_gt_atTop 0] with T hT
    have hb := norm_perron_F_horizontal_integral_le hx hσ₁ hσ₀ (neg_ne_zero.mpr hT.ne')
    rwa [neg_sq] at hb
  exact squeeze_zero_norm' h hg

#print axioms perron_horizontal_edge_bottom_tendsto_zero
