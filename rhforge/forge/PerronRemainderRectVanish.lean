import PerronRemainderDiffAt
import Mathlib.Analysis.Complex.CauchyIntegral

set_option maxHeartbeats 1000000

open Complex Filter Topology Set

open scoped Interval

/-!
# §1664 — the Perron remainder's rectangle integral VANISHES on zero-free rectangles

The Cauchy-Goursat half of the Perron contour shift: on any closed rectangle avoiding `0`, `-1`,
and the zeros of `ζ` (`s = 1` is allowed — the remainder `R = dslope Ĥ 1` is REMOVABLE there,
§1615/§1621), the boundary integral of `R` is `0`. Combined with §1662 (the pole term deposits
`2π·x²/2` crossing the same rectangle) and the §1624 split `F = pole + R`, this is exactly the
statement that shifting the Perron line across `s = 1` through a zero-free box changes `∫F` by the
residue alone — the finite-`T` contour-shift identity, modulo edge bookkeeping.
-/

theorem perron_remainder_rect_vanish {x : ℝ} (hx : 0 < x) (z w : ℂ)
    (hcond : ∀ u ∈ [[z.re, w.re]] ×ℂ [[z.im, w.im]],
      (u ≠ 1 → riemannZeta u ≠ 0) ∧ u ≠ 0 ∧ u ≠ -1) :
    (∫ t in z.re..w.re,
        dslope (Function.update
          (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
          1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1 ((t : ℂ) + z.im * Complex.I)) -
      (∫ t in z.re..w.re,
        dslope (Function.update
          (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
          1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1 ((t : ℂ) + w.im * Complex.I)) +
      Complex.I • (∫ y in z.im..w.im,
        dslope (Function.update
          (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
          1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1 ((w.re : ℂ) + y * Complex.I)) -
      Complex.I • (∫ y in z.im..w.im,
        dslope (Function.update
          (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
          1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1 ((z.re : ℂ) + y * Complex.I)) = 0 := by
  have h := integral_boundary_rect_eq_zero_of_differentiableOn
    (dslope (Function.update
      (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
      1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1) z w
    (fun u hu => (differentiableAt_perron_remainder hx (hcond u hu).1 (hcond u hu).2.1
      (hcond u hu).2.2).differentiableWithinAt)
  simpa using h

#print axioms perron_remainder_rect_vanish
