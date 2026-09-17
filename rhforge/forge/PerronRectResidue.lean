import PerronLineSplit
import PerronRemainderRectVanish
import RectangleResidue

set_option maxHeartbeats 4000000

open Complex MeasureTheory Set

open scoped Interval

/-!
# §1665 ★★★★ — THE FINITE-T PERRON SHIFT IDENTITY: `∮_∂R F = 2πi·(x²/2)` on zero-free rectangles

THE capstone of the winding arc: the full Perron integrand `F(s) = -ζ'/ζ(s)·x^{s+1}/(s(s+1))`
integrated around ANY rectangle that (a) strictly contains `s = 1` and (b) avoids `0`, `-1`, and
the zeros of `ζ`, picks up exactly the residue `2πi·x²/2`. Mechanism, fully machine-checked:
edge-by-edge split `F = pole + R` (§1614/§1624 + horizontal analogues built here), the remainder's
rectangle vanishes (§1664 Cauchy-Goursat), and the pole's rectangle winds (§1658). This is the
finite-`T` contour-shift identity: moving the truncated Perron line from `σ₀ > 1` to `σ₁ < 1`
through a zero-free box costs exactly the PNT main-term residue plus horizontal edges that vanish
as `T → ∞` (§1637). -/

theorem perron_F_rect_residue {x : ℝ} (hx : 0 < x) (z w : ℂ)
    (h1re : z.re < 1) (h2re : 1 < w.re) (h1im : z.im < 0) (h2im : 0 < w.im)
    (hcond : ∀ u ∈ [[z.re, w.re]] ×ℂ [[z.im, w.im]],
      (u ≠ 1 → riemannZeta u ≠ 0) ∧ u ≠ 0 ∧ u ≠ -1) :
    (∫ t in z.re..w.re,
        -logDeriv riemannZeta ((t : ℂ) + z.im * I) * (x : ℂ) ^ (((t : ℂ) + z.im * I) + 1) /
          (((t : ℂ) + z.im * I) * (((t : ℂ) + z.im * I) + 1))) -
      (∫ t in z.re..w.re,
        -logDeriv riemannZeta ((t : ℂ) + w.im * I) * (x : ℂ) ^ (((t : ℂ) + w.im * I) + 1) /
          (((t : ℂ) + w.im * I) * (((t : ℂ) + w.im * I) + 1))) +
      I • (∫ t in z.im..w.im,
        -logDeriv riemannZeta ((w.re : ℂ) + (t : ℝ) * I) *
            (x : ℂ) ^ (((w.re : ℂ) + (t : ℝ) * I) + 1) /
          (((w.re : ℂ) + (t : ℝ) * I) * (((w.re : ℂ) + (t : ℝ) * I) + 1))) -
      I • (∫ t in z.im..w.im,
        -logDeriv riemannZeta ((z.re : ℂ) + (t : ℝ) * I) *
            (x : ℂ) ^ (((z.re : ℂ) + (t : ℝ) * I) + 1) /
          (((z.re : ℂ) + (t : ℝ) * I) * (((z.re : ℂ) + (t : ℝ) * I) + 1)))
      = 2 * (Real.pi : ℂ) * I * ((x : ℂ) ^ (2 : ℂ) / 2) := by
  -- membership of the four edges in the closed rectangle
  have hmem_h : ∀ (c : ℝ), c ∈ [[z.im, w.im]] → ∀ t ∈ Set.uIcc z.re w.re,
      ((t : ℂ) + (c : ℝ) * I) ∈ [[z.re, w.re]] ×ℂ [[z.im, w.im]] := by
    intro c hc t ht
    rw [Complex.mem_reProdIm]
    constructor <;> simp [ht, hc]
  have hmem_v : ∀ (a : ℝ), a ∈ [[z.re, w.re]] → ∀ t ∈ Set.uIcc z.im w.im,
      ((a : ℂ) + (t : ℝ) * I) ∈ [[z.re, w.re]] ×ℂ [[z.im, w.im]] := by
    intro a ha t ht
    rw [Complex.mem_reProdIm]
    constructor <;> simp [ht, ha]
  -- horizontal split: F = pole + R along Im = c, for c ≠ 0
  have hhsplit : ∀ (c : ℝ), c ≠ 0 → c ∈ [[z.im, w.im]] →
      (∫ t in z.re..w.re,
          -logDeriv riemannZeta ((t : ℂ) + (c : ℝ) * I) *
              (x : ℂ) ^ (((t : ℂ) + (c : ℝ) * I) + 1) /
            (((t : ℂ) + (c : ℝ) * I) * (((t : ℂ) + (c : ℝ) * I) + 1)))
        = (∫ t in z.re..w.re, (x : ℂ) ^ (2 : ℂ) / 2 / (((t : ℂ) + (c : ℝ) * I) - 1))
          + (∫ t in z.re..w.re,
              dslope (Function.update
                (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
                1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1 ((t : ℂ) + (c : ℝ) * I)) := by
    intro c hc hcmem
    have hne1 : ∀ t : ℝ, ((t : ℂ) + (c : ℝ) * I) ≠ 1 := by
      intro t h
      have := congrArg Complex.im h
      simp at this
      exact hc this
    have hpole_int : IntervalIntegrable
        (fun t : ℝ => (x : ℂ) ^ (2 : ℂ) / 2 / (((t : ℂ) + (c : ℝ) * I) - 1))
        volume z.re w.re := by
      apply Continuous.intervalIntegrable
      apply Continuous.div continuous_const
      · fun_prop
      · intro t
        exact sub_ne_zero.mpr (hne1 t)
    have hR_int : IntervalIntegrable
        (fun t : ℝ =>
          dslope (Function.update
            (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
            1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1 ((t : ℂ) + (c : ℝ) * I))
        volume z.re w.re := by
      apply ContinuousOn.intervalIntegrable
      intro t ht
      have hL : ContinuousAt (fun t : ℝ => (t : ℂ) + (c : ℝ) * I) t := by fun_prop
      obtain ⟨hζ, h0, hm1⟩ := hcond ((t : ℂ) + (c : ℝ) * I) (hmem_h c hcmem t ht)
      have hR := differentiableAt_perron_remainder hx (z := (t : ℂ) + (c : ℝ) * I) hζ h0 hm1
      have hcomp : ContinuousAt
          (fun t : ℝ =>
            dslope (Function.update
              (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
              1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1 ((t : ℂ) + (c : ℝ) * I)) t :=
        ContinuousAt.comp hR.continuousAt hL
      exact hcomp.continuousWithinAt
    rw [← intervalIntegral.integral_add hpole_int hR_int]
    apply intervalIntegral.integral_congr
    intro t _ht
    have h14 := perron_dslope_eq (x := x) (hne1 t)
    simp only [Pi.add_apply]
    rw [h14]; ring
  -- the four hne packages for the two vertical §1624 applications
  have hne_r : ∀ t ∈ Set.uIcc z.im w.im,
      riemannZeta ((w.re : ℂ) + (t : ℝ) * I) ≠ 0 ∧
        ((w.re : ℂ) + (t : ℝ) * I) ≠ 0 ∧ ((w.re : ℂ) + (t : ℝ) * I) ≠ -1 := by
    intro t ht
    have hu := hcond ((w.re : ℂ) + (t : ℝ) * I) (hmem_v w.re right_mem_uIcc t ht)
    refine ⟨hu.1 ?_, hu.2.1, hu.2.2⟩
    intro h
    have h' := congrArg Complex.re h
    simp at h'
    linarith
  have hne_l : ∀ t ∈ Set.uIcc z.im w.im,
      riemannZeta ((z.re : ℂ) + (t : ℝ) * I) ≠ 0 ∧
        ((z.re : ℂ) + (t : ℝ) * I) ≠ 0 ∧ ((z.re : ℂ) + (t : ℝ) * I) ≠ -1 := by
    intro t ht
    have hu := hcond ((z.re : ℂ) + (t : ℝ) * I) (hmem_v z.re left_mem_uIcc t ht)
    refine ⟨hu.1 ?_, hu.2.1, hu.2.2⟩
    intro h
    have h' := congrArg Complex.re h
    simp at h'
    linarith
  -- split all four edges
  rw [perron_line_split hx (ne_of_gt h2re) hne_r, perron_line_split hx (ne_of_lt h1re) hne_l,
    hhsplit z.im (ne_of_lt h1im) left_mem_uIcc, hhsplit w.im (ne_of_gt h2im) right_mem_uIcc]
  -- remainder rectangle vanishes; pole rectangle winds
  have hRrect := perron_remainder_rect_vanish hx z w hcond
  have hPrect := integral_boundary_rect_residue ((x : ℂ) ^ (2 : ℂ) / 2) 1 z w
    (by simpa using h1re) (by simpa using h2re) (by simpa using h1im) (by simpa using h2im)
  simp only [smul_eq_mul, div_eq_mul_inv] at hRrect hPrect ⊢
  linear_combination hRrect + hPrect

#print axioms perron_F_rect_residue
