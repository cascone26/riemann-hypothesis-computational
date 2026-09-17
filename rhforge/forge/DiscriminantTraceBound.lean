import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Tactic

set_option maxHeartbeats 1000000

noncomputable section

/-!
# §1726 ★★★★★ — THE DISCRIMINANT HINGE: positivity of a quadratic ⟹ trace bound `t² ≤ 4q`

The second half of the function-field template's connector (companion to §1725). Weil's positivity —
the Rosati/adjunction form on correspondences is positive, which our Hodge-index negativity-on-the-
primitive-part (§1720/§1713/§1725) supplies — says a certain quadratic `n ↦ a·n² + b·n + c` in a real
parameter is `≥ 0`. THIS brick is the elementary consequence: a nonnegative quadratic with positive
leading coefficient has non-positive discriminant, `b² ≤ 4ac`. In the correspondence application
`a, b, c` are intersection numbers with `b = ±2t` (trace of Frobenius) and `ac = q·(…)`, so this
becomes exactly the trace bound `t² ≤ 4q` that §1724 turns into `|λ| = √q`.

Chain, fully abstract on our side (only the curve-numeric identification `a,b,c ↔ trace,q` is
Report-tier): Hodge negativity (§1720/§1713) → hyperbolic inequality / positivity (§1725) →
**discriminant bound `b²≤4ac`** (this brick) → trace bound `t²≤4q` → `|λ|=√q` (§1724).
-/

/-- **Discriminant hinge.** If `a·x² + b·x + c ≥ 0` for every real `x` and `a > 0`, then the
discriminant is non-positive: `b² ≤ 4·a·c`. (Evaluate at the vertex `x = −b/(2a)`.) -/
theorem disc_nonpos_of_quadratic_nonneg (a b c : ℝ) (ha : 0 < a)
    (h : ∀ x : ℝ, 0 ≤ a * x ^ 2 + b * x + c) :
    b ^ 2 ≤ 4 * a * c := by
  set x₀ : ℝ := -b / (2 * a) with hx0def
  have hne : (2 * a) ≠ 0 := ne_of_gt (by positivity)
  have hx0 : 2 * a * x₀ = -b := by rw [hx0def]; field_simp
  have hsq : (2 * a * x₀ + b) ^ 2 = 0 := by rw [hx0]; ring
  have h4a : (0:ℝ) ≤ 4 * a := by positivity
  have hpos : 0 ≤ 4 * a * (a * x₀ ^ 2 + b * x₀ + c) := mul_nonneg h4a (h x₀)
  have expand : 4 * a * (a * x₀ ^ 2 + b * x₀ + c)
      = (2 * a * x₀ + b) ^ 2 + (4 * a * c - b ^ 2) := by ring
  rw [expand, hsq] at hpos
  linarith

#print axioms disc_nonpos_of_quadratic_nonneg
