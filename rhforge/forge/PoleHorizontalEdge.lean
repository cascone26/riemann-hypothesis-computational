import PoleHorizontalBound
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic

set_option maxHeartbeats 1000000

open Complex

/-- §1635 — the pole term's horizontal-edge integral is `O(1/|T|)`:
`‖∫_{σ₁}^{σ₀} (x²/2)/((σ'+iT)-1) dσ'‖ ≤ (x²/2)/|T| · |σ₀-σ₁|`.  The per-point bound §1633 is uniform in `σ'`
(depends only on `T`), so the edge integral inherits the `1/|T|` decay — the pole term's half of the rectangle's
top/bottom edges vanishing as `|T|→∞`. -/
theorem norm_pole_horizontal_integral_le {x : ℝ} (hx : 0 < x) {σ₁ σ₀ T : ℝ} (hT : T ≠ 0) :
    ‖∫ σ' : ℝ in σ₁..σ₀, (x : ℂ) ^ (2 : ℂ) / 2 / (((σ' : ℂ) + (T : ℝ) * I) - 1)‖
      ≤ (x ^ (2 : ℝ) / 2 / |T|) * |σ₀ - σ₁| :=
  intervalIntegral.norm_integral_le_of_norm_le_const (fun σ' _ => norm_pole_le_horizontal hx hT)

#print axioms norm_pole_horizontal_integral_le
