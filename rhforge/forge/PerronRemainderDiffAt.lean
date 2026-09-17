import PerronRemainderAnalytic
import Mathlib.NumberTheory.LSeries.Dirichlet

set_option maxHeartbeats 1000000

/-!
# §1621 — pointwise holomorphy of the remainder `R = dslope Ĥ 1` (`DifferentiableAt`)

The `DifferentiableAt` (not merely `DifferentiableOn`) form of §1619: at any single `z` with `z ≠ 1 → ζ z ≠ 0`,
`z ≠ 0`, `z ≠ -1`, the remainder `R = dslope Ĥ 1` is complex-differentiable — at `z = 1` by §1615 (removable),
elsewhere by `differentiableAt_dslope_of_ne`.  This is the pointwise input to `ContinuousAt` (hence
interval-integrability of `t ↦ R(σ+it)` along vertical segments) needed to shift the Perron line integral across the
zero-free region.  EXIT0 / std-axioms.
-/

open Complex Filter Topology

theorem differentiableAt_perron_remainder {x : ℝ} (hx : 0 < x) {z : ℂ}
    (hζ : z ≠ 1 → riemannZeta z ≠ 0) (h0 : z ≠ 0) (hm1 : z ≠ -1) :
    DifferentiableAt ℂ
      (dslope
        (Function.update
          (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
          1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1) z := by
  set Ĥ : ℂ → ℂ :=
    Function.update
      (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
      1 ((x : ℂ) ^ (2 : ℂ) / 2) with hĤ
  rcases eq_or_ne z 1 with rfl | hz1
  · exact (analyticAt_perron_remainder hx).differentiableAt
  · have heq : Ĥ =ᶠ[nhds z]
        (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1)))) := by
      filter_upwards [isOpen_ne.mem_nhds hz1] with u hu; rw [hĤ]; exact Function.update_of_ne hu _ _
    have hĤd : DifferentiableAt ℂ Ĥ z := by
      rw [heq.differentiableAt_iff]
      exact (differentiableAt_id.sub_const 1).mul
        (differentiableAt_perron_integrand hx hz1 (hζ hz1) h0 hm1)
    exact (differentiableAt_dslope_of_ne hz1).mpr hĤd

#print axioms differentiableAt_perron_remainder
