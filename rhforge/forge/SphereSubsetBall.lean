import Mathlib.Analysis.SpecialFunctions.Complex.Circle

/-!
# §1397 — `sphere c R ⊆ ball c' R''` when `dist c c' + R < R''`

Off-center sphere-in-ball inclusion (triangle inequality).  Lets a ζ-bound proven on a ball centred at `1+iγ`
(§1130) transfer to the sphere centred at the test point `σ+iγ`: with `dist(σ+iγ, 1+iγ) = σ-1` and `R'' > σ-1+R`,
`sphere (σ+iγ) R ⊆ ball (1+iγ) R''`, so `‖ζ‖ ≤ Mζ` on the sphere.  Supplies §1378's `hζsphere`.  EXIT0 / std.
-/

open Metric

theorem sphere_subset_ball' {c c' : ℂ} {R R'' : ℝ} (h : dist c c' + R < R'') :
    sphere c R ⊆ ball c' R'' := by
  intro s hs
  rw [mem_sphere] at hs
  rw [mem_ball]
  calc dist s c' ≤ dist s c + dist c c' := dist_triangle _ _ _
    _ = R + dist c c' := by rw [hs]
    _ < R'' := by linarith

#print axioms sphere_subset_ball'
