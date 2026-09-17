import PerronRectResidue

set_option maxHeartbeats 4000000

open Complex MeasureTheory Set

open scoped Interval

/-!
# §1666 — the Perron VERTICAL-LINE SHIFT form: `∫_{Re=σ₀} F - ∫_{Re=σ₁} F = 2π·x²/2 + connectors`

§1665 rearranged into contour-shift form (the same move §1660 made for the bare pole): the
difference of the two truncated vertical Perron line integrals across `s = 1` equals the clean
main-term deposit `2π·(x²/2)` plus the two horizontal connector integrals of `F` — which vanish
as `T → ∞` by §1630/§1637's uniform bound. The direct input for the effective-PNT assembly:
pick the zero-free box from the DLVP region, shift, bound the connectors.
-/

theorem perron_F_vline_shift {x : ℝ} (hx : 0 < x) (z w : ℂ)
    (h1re : z.re < 1) (h2re : 1 < w.re) (h1im : z.im < 0) (h2im : 0 < w.im)
    (hcond : ∀ u ∈ [[z.re, w.re]] ×ℂ [[z.im, w.im]],
      (u ≠ 1 → riemannZeta u ≠ 0) ∧ u ≠ 0 ∧ u ≠ -1) :
    (∫ t in z.im..w.im,
        -logDeriv riemannZeta ((w.re : ℂ) + (t : ℝ) * I) *
            (x : ℂ) ^ (((w.re : ℂ) + (t : ℝ) * I) + 1) /
          (((w.re : ℂ) + (t : ℝ) * I) * (((w.re : ℂ) + (t : ℝ) * I) + 1))) -
      (∫ t in z.im..w.im,
        -logDeriv riemannZeta ((z.re : ℂ) + (t : ℝ) * I) *
            (x : ℂ) ^ (((z.re : ℂ) + (t : ℝ) * I) + 1) /
          (((z.re : ℂ) + (t : ℝ) * I) * (((z.re : ℂ) + (t : ℝ) * I) + 1)))
      = 2 * (Real.pi : ℂ) * ((x : ℂ) ^ (2 : ℂ) / 2) +
        I * ((∫ t in z.re..w.re,
            -logDeriv riemannZeta ((t : ℂ) + z.im * I) * (x : ℂ) ^ (((t : ℂ) + z.im * I) + 1) /
              (((t : ℂ) + z.im * I) * (((t : ℂ) + z.im * I) + 1))) -
          (∫ t in z.re..w.re,
            -logDeriv riemannZeta ((t : ℂ) + w.im * I) * (x : ℂ) ^ (((t : ℂ) + w.im * I) + 1) /
              (((t : ℂ) + w.im * I) * (((t : ℂ) + w.im * I) + 1)))) := by
  have h := perron_F_rect_residue hx z w h1re h2re h1im h2im hcond
  simp only [smul_eq_mul] at h
  linear_combination (-I) * h +
    ((∫ t in z.im..w.im,
        -logDeriv riemannZeta ((w.re : ℂ) + (t : ℝ) * I) *
            (x : ℂ) ^ (((w.re : ℂ) + (t : ℝ) * I) + 1) /
          (((w.re : ℂ) + (t : ℝ) * I) * (((w.re : ℂ) + (t : ℝ) * I) + 1))) -
      (∫ t in z.im..w.im,
        -logDeriv riemannZeta ((z.re : ℂ) + (t : ℝ) * I) *
            (x : ℂ) ^ (((z.re : ℂ) + (t : ℝ) * I) + 1) /
          (((z.re : ℂ) + (t : ℝ) * I) * (((z.re : ℂ) + (t : ℝ) * I) + 1))) -
      2 * (Real.pi : ℂ) * ((x : ℂ) ^ (2 : ℂ) / 2)) * Complex.I_sq

#print axioms perron_F_vline_shift
