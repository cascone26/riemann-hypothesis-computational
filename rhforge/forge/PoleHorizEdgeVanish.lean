import PoleHorizontalEdge
import Mathlib.Topology.Order.Bornology

set_option maxHeartbeats 1000000

open Complex Filter Topology

/-- §1638 — the pole term's horizontal edge also VANISHES as `T → ∞`:
`∫_{σ₁}^{σ₀} (x²/2)/((σ'+iT)-1) dσ' → 0`.  Squeeze the §1635 bound `≤ (x²/2)/|T|·|σ₀-σ₁| = O(1/|T|)`.
Companion of §1637; together the whole integrand's top/bottom rectangle edges disappear at infinity. -/
theorem pole_horizontal_edge_tendsto_zero {x σ₁ σ₀ : ℝ} (hx : 0 < x) :
    Tendsto (fun T : ℝ => ∫ σ' : ℝ in σ₁..σ₀,
        (x : ℂ) ^ (2 : ℂ) / 2 / (((σ' : ℂ) + (T : ℝ) * I) - 1)) atTop (𝓝 0) := by
  have hg : Tendsto (fun T : ℝ => x ^ (2 : ℝ) / 2 / |T| * |σ₀ - σ₁|) atTop (𝓝 0) := by
    have h1 : Tendsto (fun T : ℝ => x ^ (2 : ℝ) / 2 / |T|) atTop (𝓝 0) :=
      Tendsto.div_atTop tendsto_const_nhds tendsto_abs_atTop_atTop
    simpa using h1.mul_const |σ₀ - σ₁|
  have h : ∀ᶠ T : ℝ in atTop, ‖∫ σ' : ℝ in σ₁..σ₀,
        (x : ℂ) ^ (2 : ℂ) / 2 / (((σ' : ℂ) + (T : ℝ) * I) - 1)‖
      ≤ x ^ (2 : ℝ) / 2 / |T| * |σ₀ - σ₁| := by
    filter_upwards [eventually_ne_atTop 0] with T hT
    exact norm_pole_horizontal_integral_le hx hT
  exact squeeze_zero_norm' h hg

#print axioms pole_horizontal_edge_tendsto_zero
