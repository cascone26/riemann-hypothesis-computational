import WindingEdges
import ArctanCrossSum

set_option maxHeartbeats 4000000

open Complex intervalIntegral Real

/-!
# §1657 ★★★ — THE RECTANGLE-WINDING THEOREM: `∮_∂R (z - p)⁻¹ dz = 2πi` for `p` inside `R`

THE lemma Mathlib LACKS (it has circle winding `circleIntegral.integral_sub_inv_of_mem_ball` and
rectangle Cauchy-Goursat VANISHING only). Stated in Mathlib's own rectangle-boundary form (the exact
4-integral combination of `integral_boundary_rect_eq_zero_of_differentiableOn`), so it composes
directly with Cauchy-Goursat to give the Cauchy integral formula / residue deposit on rectangles.
Proof: the four §1655 edge closed forms; the eight ½log corner terms cancel in pairs (adjacent edges
meet at corners), and the eight arctan terms are the §1656 cross sum = 2π. Branch-cut-free.
-/

theorem integral_boundary_rect_sub_inv (p z w : ℂ)
    (h1 : z.re < p.re) (h2 : p.re < w.re) (h3 : z.im < p.im) (h4 : p.im < w.im) :
    (∫ x in z.re..w.re, ((x : ℂ) + z.im * Complex.I - p)⁻¹) -
      (∫ x in z.re..w.re, ((x : ℂ) + w.im * Complex.I - p)⁻¹) +
      Complex.I • (∫ y in z.im..w.im, ((w.re : ℂ) + y * Complex.I - p)⁻¹) -
      Complex.I • (∫ y in z.im..w.im, ((z.re : ℂ) + y * Complex.I - p)⁻¹)
      = 2 * (Real.pi : ℂ) * Complex.I := by
  have hbot : z.im ≠ p.im := ne_of_lt h3
  have htop : w.im ≠ p.im := ne_of_gt h4
  have hright : w.re ≠ p.re := ne_of_gt h2
  have hleft : z.re ≠ p.re := ne_of_lt h1
  have hA : 0 < w.re - p.re := sub_pos.mpr h2
  have hB : 0 < p.re - z.re := sub_pos.mpr h1
  have hC : 0 < w.im - p.im := sub_pos.mpr h4
  have hD : 0 < p.im - z.im := sub_pos.mpr h3
  -- the four edge closed forms
  rw [integral_horizontal_winding p hbot, integral_horizontal_winding p htop,
    integral_vertical_winding p hright, integral_vertical_winding p hleft]
  simp only [smul_eq_mul]
  unfold windingPrim
  -- normalize the four vertical-edge log arguments to the horizontal-edge corner form
  have l1 : (z.im - p.im) ^ 2 + (p.re - w.re) ^ 2 = (w.re - p.re) ^ 2 + (z.im - p.im) ^ 2 := by
    ring
  have l2 : (z.im - p.im) ^ 2 + (p.re - z.re) ^ 2 = (z.re - p.re) ^ 2 + (z.im - p.im) ^ 2 := by
    ring
  have l3 : (w.im - p.im) ^ 2 + (p.re - w.re) ^ 2 = (w.re - p.re) ^ 2 + (w.im - p.im) ^ 2 := by
    ring
  have l4 : (w.im - p.im) ^ 2 + (p.re - z.re) ^ 2 = (z.re - p.re) ^ 2 + (w.im - p.im) ^ 2 := by
    ring
  rw [l1, l2, l3, l4]
  -- normalize the eight arctan arguments to the canonical positive ratios A,B,C,D
  have rb1 : (w.re - p.re) / (z.im - p.im) = -((w.re - p.re) / (p.im - z.im)) := by
    rw [show z.im - p.im = -(p.im - z.im) by ring, div_neg]
  have rb2 : (z.re - p.re) / (z.im - p.im) = (p.re - z.re) / (p.im - z.im) := by
    rw [show z.re - p.re = -(p.re - z.re) by ring, show z.im - p.im = -(p.im - z.im) by ring,
      neg_div_neg_eq]
  have rt2 : (z.re - p.re) / (w.im - p.im) = -((p.re - z.re) / (w.im - p.im)) := by
    rw [show z.re - p.re = -(p.re - z.re) by ring, neg_div]
  have rr1 : (w.im - p.im) / (p.re - w.re) = -((w.im - p.im) / (w.re - p.re)) := by
    rw [show p.re - w.re = -(w.re - p.re) by ring, div_neg]
  have rr2 : (z.im - p.im) / (p.re - w.re) = (p.im - z.im) / (w.re - p.re) := by
    rw [show z.im - p.im = -(p.im - z.im) by ring, show p.re - w.re = -(w.re - p.re) by ring,
      neg_div_neg_eq]
  have rl2 : (z.im - p.im) / (p.re - z.re) = -((p.im - z.im) / (p.re - z.re)) := by
    rw [show z.im - p.im = -(p.im - z.im) by ring, neg_div]
  rw [rb1, rb2, rt2, rr1, rr2, rl2]
  simp only [Real.arctan_neg]
  -- the winding angle: §1656 cross sum, cast to ℂ
  have hsum := arctan_rectangle_cross_sum hA hB hC hD
  have hsumC := congrArg (fun r : ℝ => (r : ℂ)) hsum
  push_cast at hsumC
  push_cast
  linear_combination Complex.I * hsumC +
    ((Real.log ((w.re - p.re) ^ 2 + (z.im - p.im) ^ 2) : ℂ) * (1 / 2) -
      (Real.log ((z.re - p.re) ^ 2 + (z.im - p.im) ^ 2) : ℂ) * (1 / 2) -
      (Real.log ((w.re - p.re) ^ 2 + (w.im - p.im) ^ 2) : ℂ) * (1 / 2) +
      (Real.log ((z.re - p.re) ^ 2 + (w.im - p.im) ^ 2) : ℂ) * (1 / 2) -
      Complex.I * ((Real.arctan ((w.im - p.im) / (w.re - p.re)) : ℂ) +
        (Real.arctan ((p.im - z.im) / (w.re - p.re)) : ℂ) +
        (Real.arctan ((w.im - p.im) / (p.re - z.re)) : ℂ) +
        (Real.arctan ((p.im - z.im) / (p.re - z.re)) : ℂ))) * Complex.I_sq

#print axioms integral_boundary_rect_sub_inv
