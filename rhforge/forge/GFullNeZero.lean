import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Analysis.Meromorphic.Divisor

/-!
# §1382 — `g_full = P_rest · g` is nonvanishing on the inner disk

In the two-radii split, `P_rest = ∏_{ρ∈T} (s-ρ)^{m ρ}` carries the *annulus* zeros (`ρ ∉ closedBall c (R/2)`).
On the inner disk `ball c (R/2)` every such factor is nonzero — a point `s` with `‖s-c‖ < R/2` cannot equal a
zero `ρ` with `‖ρ-c‖ > R/2` — and `g` is nonvanishing everywhere (§1343), so `g_full = P_rest · g` is nonvanishing
on `ball c (R/2)`.  This is the `hgnz` hypothesis §1378 needs.  EXIT0 / std-axioms.
-/

theorem prod_rest_mul_ne_zero {T : Finset ℂ} {c s : ℂ} {R : ℝ} (m : ℂ → ℕ) (g : ℂ → ℂ)
    (hT : ∀ ρ ∈ T, ρ ∉ Metric.closedBall c (R / 2))
    (hs : s ∈ Metric.ball c (R / 2)) (hg : g s ≠ 0) :
    (∏ ρ ∈ T, (s - ρ) ^ (m ρ)) * g s ≠ 0 := by
  apply mul_ne_zero _ hg
  apply Finset.prod_ne_zero_iff.mpr
  intro ρ hρ
  apply pow_ne_zero
  rw [sub_ne_zero]
  intro h
  exact hT ρ hρ (h ▸ Metric.ball_subset_closedBall hs)

#print axioms prod_rest_mul_ne_zero
