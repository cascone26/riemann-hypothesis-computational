import Mathlib.Analysis.SpecialFunctions.Trigonometric.Arctan
import Mathlib.Tactic

set_option maxHeartbeats 1000000

open Real

/-!
# §1656 — the ARCTAN CROSS SUM: the winding angle of a rectangle is `2π`

The angular content of the rectangle-winding lemma, isolated as pure real analysis: at each of the
four corners of a rectangle strictly containing the pole, the two adjacent edges contribute a pair
`arctan(u/v) + arctan(v/u) = π/2` (one aspect ratio from the horizontal edge, its reciprocal from
the vertical edge). Four corners × `π/2` = `2π` — the total winding angle.
-/

/-- The corner pair: `arctan(x/y) + arctan(y/x) = π/2` for `x, y > 0`
(reciprocal-argument complementarity, from `arctan_inv_of_pos`). -/
theorem arctan_div_add_arctan_div_swap {x y : ℝ} (hx : 0 < x) (hy : 0 < y) :
    arctan (x / y) + arctan (y / x) = π / 2 := by
  have h := arctan_inv_of_pos (div_pos hx hy)
  rw [inv_div] at h
  linarith

/-- §1656 — the full rectangle winding angle: with the pole at distances `A` (right), `B` (left),
`C` (up), `D` (down) from the four edges (all positive = pole strictly inside), the eight arctan
terms produced by the four `windingPrim` edge evaluations sum to `2π`. -/
theorem arctan_rectangle_cross_sum {A B C D : ℝ}
    (hA : 0 < A) (hB : 0 < B) (hC : 0 < C) (hD : 0 < D) :
    (arctan (A / C) + arctan (C / A)) + (arctan (B / C) + arctan (C / B))
      + (arctan (A / D) + arctan (D / A)) + (arctan (B / D) + arctan (D / B)) = 2 * π := by
  have h1 := arctan_div_add_arctan_div_swap hA hC
  have h2 := arctan_div_add_arctan_div_swap hB hC
  have h3 := arctan_div_add_arctan_div_swap hA hD
  have h4 := arctan_div_add_arctan_div_swap hB hD
  linarith

#print axioms arctan_div_add_arctan_div_swap
#print axioms arctan_rectangle_cross_sum
