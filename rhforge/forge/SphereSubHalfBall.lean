import Mathlib.Analysis.SpecialFunctions.Complex.Circle

/-!
# §1374 — half-disk zeros stay `≥ R/2` from the `R`-sphere

The geometric key to the factor sup-bound: a zero `ρ` inside the *half* disk `closedBall c (R/2)` is at distance
`≥ R/2` from every point `s` of the boundary sphere `sphere c R`, since
`‖s - ρ‖ ≥ ‖s - c‖ - ‖ρ - c‖ = R - ‖ρ - c‖ ≥ R - R/2 = R/2`.

So the partial product over half-disk zeros, `∏ (s - ρ)`, is `≥ (R/2)^n` on the sphere (`n` = zero count,
bounded by Jensen §1363) — the `Plow` that §1373 turns into `‖g‖ ≤ Mζ/Plow`.  This is the standard two-radii
construction of the quantitative zero-free region.  EXIT0 / std-axioms.
-/

open Complex

theorem sphere_sub_halfball_ge {c s ρ : ℂ} {R : ℝ}
    (hs : s ∈ Metric.sphere c R) (hρ : ρ ∈ Metric.closedBall c (R / 2)) :
    R / 2 ≤ ‖s - ρ‖ := by
  rw [Metric.mem_sphere, Complex.dist_eq] at hs
  rw [Metric.mem_closedBall, Complex.dist_eq] at hρ
  have h := norm_sub_norm_le (s - c) (ρ - c)
  have heq : s - c - (ρ - c) = s - ρ := by ring
  rw [heq] at h
  linarith

#print axioms sphere_sub_halfball_ge
