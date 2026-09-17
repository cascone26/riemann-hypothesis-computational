import PerronRemainderDiffOn
import Mathlib.Analysis.Complex.CauchyIntegral

set_option maxHeartbeats 1000000

/-!
# §1620 — the remainder's rectangle boundary integral vanishes (even across `s = 1`)

Since the remainder `R = dslope Ĥ 1` is holomorphic on the whole closed rectangle whenever that rectangle is free of
`ζ`-zeros (except possibly `1`, where `R` is regular) and avoids `{0,-1}` (§1619), Cauchy–Goursat
(`integral_boundary_rect_eq_zero_of_differentiableOn`) makes its boundary integral `0` — unlike the raw integrand
`F`, this holds even for a rectangle straddling `s=1`.  Combined with the split §1614 (`F = (x²/2)/(s-1) + R`), the
entire contour contribution of `F` reduces to the explicit pole term, whose value is the §1618 residue `2πi·x²/2`.
EXIT0 / std-axioms.
-/

open Complex

theorem perron_remainder_rect_boundary_zero {x : ℝ} (hx : 0 < x) {z w : ℂ}
    (hζ : ∀ p ∈ (Complex.reProdIm (Set.uIcc z.re w.re) (Set.uIcc z.im w.im)),
      p ≠ 1 → riemannZeta p ≠ 0)
    (h0 : ∀ p ∈ (Complex.reProdIm (Set.uIcc z.re w.re) (Set.uIcc z.im w.im)), p ≠ 0)
    (hm1 : ∀ p ∈ (Complex.reProdIm (Set.uIcc z.re w.re) (Set.uIcc z.im w.im)), p ≠ -1) :
    (∫ t : ℝ in z.re..w.re,
        dslope (Function.update
          (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
          1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1 (↑t + ↑z.im * I)) -
      (∫ t : ℝ in z.re..w.re,
        dslope (Function.update
          (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
          1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1 (↑t + ↑w.im * I)) +
      I • (∫ y : ℝ in z.im..w.im,
        dslope (Function.update
          (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
          1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1 (↑w.re + ↑y * I)) -
      I • (∫ y : ℝ in z.im..w.im,
        dslope (Function.update
          (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
          1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1 (↑z.re + ↑y * I)) = 0 := by
  refine integral_boundary_rect_eq_zero_of_differentiableOn
    (dslope (Function.update
      (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
      1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1) z w ?_
  exact differentiableOn_perron_remainder hx hζ h0 hm1

#print axioms perron_remainder_rect_boundary_zero
