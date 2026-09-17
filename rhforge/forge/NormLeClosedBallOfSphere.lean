import Mathlib.Analysis.Complex.AbsMax

/-!
# §1370 — max-modulus on a ball: `sup` over the closed disk `≤ sup` over the sphere

The maximum-modulus principle in the shape the factor bound needs: if `g` is holomorphic on `ball c R` and
continuous up to the boundary, and `‖g‖ ≤ C` on the bounding sphere, then `‖g‖ ≤ C` throughout `closedBall c R`.
A direct specialisation of Mathlib's `Complex.norm_le_of_forall_mem_frontier_norm_le` (frontier of the ball is the
sphere, closure is the closed ball).  Applied to the nonvanishing factor `g = ζ/P`, this reduces the disk bound
`sup‖ζ/P‖ ≤ M'` (needed by §1369 for `hgB`) to a sphere bound.  EXIT0 / std-axioms.
-/

open Complex Metric

theorem norm_le_on_closedBall_of_sphere {c : ℂ} {R C : ℝ} (hR : 0 < R) {g : ℂ → ℂ}
    (hg : DiffContOnCl ℂ g (ball c R)) (hC : ∀ z ∈ sphere c R, ‖g z‖ ≤ C)
    {z : ℂ} (hz : z ∈ closedBall c R) : ‖g z‖ ≤ C := by
  have hCf : ∀ w ∈ frontier (ball c R), ‖g w‖ ≤ C := by
    intro w hw; rw [frontier_ball c hR.ne'] at hw; exact hC w hw
  have hzc : z ∈ closure (ball c R) := by rw [closure_ball c hR.ne']; exact hz
  exact Complex.norm_le_of_forall_mem_frontier_norm_le Metric.isBounded_ball hg hCf hzc

#print axioms norm_le_on_closedBall_of_sphere
