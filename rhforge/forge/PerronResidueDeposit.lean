import PerronRemainderAnalytic
import Mathlib.Analysis.Complex.CauchyIntegral

set_option maxHeartbeats 1000000

/-!
# §1616 — the residue deposit: `∮_{C(1,R)} Ĥ(z)/(z-1) dz = 2πi·(x²/2)`

Cauchy's integral formula for `Ĥ = update ((s-1)F) 1 (x²/2)` (analytic on a small disc around `1` by §1613) gives
`∮_{C(1,R)} Ĥ(z)/(z-1) dz = 2πi·Ĥ(1) = 2πi·(x²/2)`.  Since `F(z) = Ĥ(z)/(z-1)` off `z=1` (and the circle avoids
`1`), this is the residue the Perron integrand deposits at its pole — the `x²/2` PNT main term, now realized as a
genuine contour integral.  The disc's freedom from `ζ`-zeros / `{0,-1}` is the geometric input (supplied downstream
from the zero-free region §1575); here it is a clean hypothesis.  EXIT0 / std-axioms.
-/

open Complex Metric Topology

theorem perron_residue_deposit {x : ℝ} (hx : 0 < x) {R : ℝ} (hR : 0 < R)
    (hζ : ∀ z ∈ closedBall (1 : ℂ) R, z ≠ 1 → riemannZeta z ≠ 0)
    (h0 : ∀ z ∈ closedBall (1 : ℂ) R, z ≠ 0)
    (hm1 : ∀ z ∈ closedBall (1 : ℂ) R, z ≠ -1) :
    (∮ z in C((1 : ℂ), R),
        (Function.update
          (fun w => (w - 1) * (-logDeriv riemannZeta w * (x : ℂ) ^ (w + 1) / (w * (w + 1))))
          1 ((x : ℂ) ^ (2 : ℂ) / 2)) z / (z - 1))
      = 2 * ↑Real.pi * I * ((x : ℂ) ^ (2 : ℂ) / 2) := by
  set Ĥ : ℂ → ℂ :=
    Function.update
      (fun w => (w - 1) * (-logDeriv riemannZeta w * (x : ℂ) ^ (w + 1) / (w * (w + 1))))
      1 ((x : ℂ) ^ (2 : ℂ) / 2) with hĤ
  -- `Ĥ` is differentiable on the whole closed disc
  have hdiff : ∀ z ∈ closedBall (1 : ℂ) R, DifferentiableAt ℂ Ĥ z := by
    intro z hz
    rcases eq_or_ne z 1 with rfl | hz1
    · exact (analyticAt_sub_one_mul_perron_integrand hx).differentiableAt
    · have heq : Ĥ =ᶠ[𝓝 z]
          (fun w => (w - 1) * (-logDeriv riemannZeta w * (x : ℂ) ^ (w + 1) / (w * (w + 1)))) := by
        filter_upwards [isOpen_ne.mem_nhds hz1] with w hw; rw [hĤ]; exact Function.update_of_ne hw _ _
      rw [heq.differentiableAt_iff]
      exact (differentiableAt_id.sub_const 1).mul
        (differentiableAt_perron_integrand hx hz1 (hζ z hz hz1) (h0 z hz) (hm1 z hz))
  have hdon : DifferentiableOn ℂ Ĥ (closedBall (1 : ℂ) R) :=
    fun z hz => (hdiff z hz).differentiableWithinAt
  have hcont : ContinuousOn Ĥ (closedBall (1 : ℂ) R) := hdon.continuousOn
  have hres := circleIntegral_div_sub_of_differentiable_on_off_countable
    (Set.countable_empty) (mem_ball_self hR) hcont
    (fun z hz => hdiff z (ball_subset_closedBall hz.1))
  rw [hres, hĤ, Function.update_self]

#print axioms perron_residue_deposit
