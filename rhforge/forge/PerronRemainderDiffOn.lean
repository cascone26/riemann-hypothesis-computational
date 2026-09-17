import PerronRemainderAnalytic
import Mathlib.NumberTheory.LSeries.Dirichlet

set_option maxHeartbeats 1000000

/-!
# §1619 — the remainder `R = dslope Ĥ 1` is holomorphic on any zero-free region (incl. across `s=1`)

On any set `U` where each point `z` satisfies `z ≠ 1 → ζ z ≠ 0`, `z ≠ 0`, `z ≠ -1`, the holomorphic remainder
`R = dslope Ĥ 1` is `DifferentiableOn`: at `z = 1` by §1615 (removable, `R` analytic there), and at `z ≠ 1` by
`differentiableAt_dslope_of_ne` reducing to `Ĥ = (z-1)F` differentiable (§1607).  Crucially this holds even when `U`
straddles `s = 1` — `R` has NO pole there — so `R`'s rectangle contour integral vanishes across the former pole,
isolating the entire residue into the explicit `(x²/2)/(s-1)` term.  EXIT0 / std-axioms.
-/

open Complex

theorem differentiableOn_perron_remainder {x : ℝ} (hx : 0 < x) {U : Set ℂ}
    (hζ : ∀ z ∈ U, z ≠ 1 → riemannZeta z ≠ 0) (h0 : ∀ z ∈ U, z ≠ 0)
    (hm1 : ∀ z ∈ U, z ≠ -1) :
    DifferentiableOn ℂ
      (dslope
        (Function.update
          (fun z => (z - 1) * (-logDeriv riemannZeta z * (x : ℂ) ^ (z + 1) / (z * (z + 1))))
          1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1) U := by
  set Ĥ : ℂ → ℂ :=
    Function.update
      (fun z => (z - 1) * (-logDeriv riemannZeta z * (x : ℂ) ^ (z + 1) / (z * (z + 1))))
      1 ((x : ℂ) ^ (2 : ℂ) / 2) with hĤ
  intro z hz
  rcases eq_or_ne z 1 with rfl | hz1
  · exact (analyticAt_perron_remainder hx).differentiableAt.differentiableWithinAt
  · have heq : Ĥ =ᶠ[nhds z]
        (fun w => (w - 1) * (-logDeriv riemannZeta w * (x : ℂ) ^ (w + 1) / (w * (w + 1)))) := by
      filter_upwards [isOpen_ne.mem_nhds hz1] with w hw; rw [hĤ]; exact Function.update_of_ne hw _ _
    have hĤd : DifferentiableAt ℂ Ĥ z := by
      rw [heq.differentiableAt_iff]
      exact (differentiableAt_id.sub_const 1).mul
        (differentiableAt_perron_integrand hx hz1 (hζ z hz hz1) (h0 z hz) (hm1 z hz))
    exact ((differentiableAt_dslope_of_ne hz1).mpr hĤd).differentiableWithinAt

#print axioms differentiableOn_perron_remainder
