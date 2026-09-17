import PerronIntegrandDiffOn
import Mathlib.Analysis.Complex.CauchyIntegral

set_option maxHeartbeats 1000000

/-!
# §1611 — Cauchy–Goursat for the Perron integrand: rectangle boundary integral vanishes (Re > 1)

For a rectangle with both vertical sides at `Re > 1`, the integrand `F(s) = -ζ'/ζ(s)·x^{s+1}/(s(s+1))` is holomorphic
on the whole closed rectangle (§1610 — no pole, no ζ-zeros there), so by Mathlib's rectangle Cauchy–Goursat
(`integral_boundary_rect_eq_zero_of_differentiableOn`) its boundary integral is `0`.  Equivalently: the smoothed
Perron integral is independent of the vertical line `σ > 1` (up to the horizontal side contributions).  This is the
first genuine contour-shift identity; crossing `s=1` to pick up the §1608 residue `x²/2` is the next (deeper) step.
EXIT0 / std-axioms.
-/

open Complex

theorem perron_rect_boundary_zero {x : ℝ} (hx : 0 < x) {z w : ℂ}
    (hz : 1 < z.re) (hw : 1 < w.re) :
    (∫ t : ℝ in z.re..w.re,
        (fun s => -logDeriv riemannZeta s * (x : ℂ) ^ (s + 1) / (s * (s + 1)))
          (↑t + ↑z.im * I)) -
      (∫ t : ℝ in z.re..w.re,
        (fun s => -logDeriv riemannZeta s * (x : ℂ) ^ (s + 1) / (s * (s + 1)))
          (↑t + ↑w.im * I)) +
      I • (∫ y : ℝ in z.im..w.im,
        (fun s => -logDeriv riemannZeta s * (x : ℂ) ^ (s + 1) / (s * (s + 1)))
          (↑w.re + ↑y * I)) -
      I • (∫ y : ℝ in z.im..w.im,
        (fun s => -logDeriv riemannZeta s * (x : ℂ) ^ (s + 1) / (s * (s + 1)))
          (↑z.re + ↑y * I)) = 0 := by
  refine integral_boundary_rect_eq_zero_of_differentiableOn
    (fun s => -logDeriv riemannZeta s * (x : ℂ) ^ (s + 1) / (s * (s + 1))) z w ?_
  apply differentiableOn_perron_integrand_reGT hx
  intro p hp
  have hp1 : (z.re ⊓ w.re) ≤ p.re := ((Complex.mem_reProdIm.mp hp).1).1
  have hlt : (1 : ℝ) < z.re ⊓ w.re := lt_inf_iff.mpr ⟨hz, hw⟩
  linarith

#print axioms perron_rect_boundary_zero
