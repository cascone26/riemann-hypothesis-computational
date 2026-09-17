import PerronIntegrandDiff
import Mathlib.NumberTheory.LSeries.Dirichlet

set_option maxHeartbeats 1000000

/-!
# §1610 — the Perron integrand is holomorphic on any region right of `Re = 1`

On any set `U` where every point has `Re > 1`, the integrand `F(s) = -ζ'/ζ(s)·x^{s+1}/(s(s+1))` is `DifferentiableOn`:
`Re s > 1` forces `s ≠ 1, 0, -1` and `ζ s ≠ 0` (`riemannZeta_ne_zero_of_one_lt_re`), so §1607 applies pointwise.
This is the direct input to the Cauchy–Goursat rectangle lemma
`integral_boundary_rect_eq_zero_of_differentiableOn` for a rectangle in the pole-free half-plane — giving contour
independence of the smoothed Perron integral between two vertical lines `σ > 1`, the first real contour-shift step.
EXIT0 / std-axioms.
-/

open Complex

theorem differentiableOn_perron_integrand_reGT {x : ℝ} (hx : 0 < x) {U : Set ℂ}
    (hU : ∀ p ∈ U, 1 < p.re) :
    DifferentiableOn ℂ
      (fun s => -logDeriv riemannZeta s * (x : ℂ) ^ (s + 1) / (s * (s + 1))) U := by
  intro p hp
  have hpre : 1 < p.re := hU p hp
  have hp1 : p ≠ 1 := by intro h; rw [h] at hpre; simp at hpre
  have hpζ : riemannZeta p ≠ 0 := riemannZeta_ne_zero_of_one_lt_re hpre
  have hp0 : p ≠ 0 := by intro h; rw [h] at hpre; norm_num at hpre
  have hpm1 : p ≠ -1 := by intro h; rw [h] at hpre; norm_num at hpre
  exact (differentiableAt_perron_integrand hx hp1 hpζ hp0 hpm1).differentiableWithinAt

#print axioms differentiableOn_perron_integrand_reGT
