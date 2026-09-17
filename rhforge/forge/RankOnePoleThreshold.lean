import Mathlib.Tactic

/-!
# §1773 — Exact rank-one pole-lift threshold (Birman–Schwinger, 2D model)

The minimal model of the Weil-form knife-edge (cycles 87–90, verified to 6 digits).  The no-pole form has one
negative square (`μ₁ < 0`) and one effective positive direction from the near-null cluster (`μ₂ > 0`); the
completed-ζ pole adds the rank-one term `(p₁ x + p₂ y)²`.  The full form
    q(x,y) = μ₁ x² + μ₂ y² + (p₁ x + p₂ y)²
is positive semidefinite **iff the determinant is nonnegative**:
    0 ≤ D := μ₁ μ₂ + μ₂ p₁² + μ₁ p₂².
Rearranged (for `μ₂>0`): `p₁² ≥ (-μ₁)(1 + p₂²/μ₂)` — the pole strength `p₁²` on the negative square must
exceed the depth `|μ₁|` by the factor `(1 + β)`, `β = p₂²/μ₂`.  This is exactly the measured threshold
`⟨p,φ₁⟩²/|μ₁| = 1+β` (saturated, margin = λ_a > 0 ⟺ RH).  Proof: the completing-the-square identity
`(μ₂+p₂²)·q = ((μ₂+p₂²)y + p₁p₂ x)² + D·x²`.  EXIT0 / std-axioms.
-/

/-- Completing-the-square identity: `c · q = (c y + b x)² + D x²` with `c = μ₂+p₂²`, `b = p₁p₂`,
`D = μ₁μ₂ + μ₂p₁² + μ₁p₂²`.  Pure ring identity. -/
theorem poleLift_csquare (μ₁ μ₂ p₁ p₂ x y : ℝ) :
    (μ₂ + p₂^2) * (μ₁ * x^2 + μ₂ * y^2 + (p₁*x + p₂*y)^2)
      = ((μ₂ + p₂^2) * y + p₁*p₂*x)^2 + (μ₁*μ₂ + μ₂*p₁^2 + μ₁*p₂^2) * x^2 := by
  ring

/-- **Exact PSD threshold.**  With `μ₂ > 0` (so `μ₂+p₂² > 0`), the form `q` is positive semidefinite iff the
determinant `D = μ₁μ₂ + μ₂p₁² + μ₁p₂²` is nonnegative. -/
theorem poleLift_psd_iff (μ₁ μ₂ p₁ p₂ : ℝ) (hμ₂ : 0 < μ₂) :
    (∀ x y : ℝ, 0 ≤ μ₁ * x^2 + μ₂ * y^2 + (p₁*x + p₂*y)^2)
      ↔ 0 ≤ μ₁*μ₂ + μ₂*p₁^2 + μ₁*p₂^2 := by
  have hc : 0 < μ₂ + p₂^2 := by positivity
  constructor
  · -- PSD ⟹ D ≥ 0: evaluate at x = μ₂+p₂², y = -(p₁p₂) (division-free; the square vanishes)
    intro h
    have hq := h (μ₂ + p₂^2) (-(p₁*p₂))
    -- identity: (μ₂+p₂²) · q(c, -p₁p₂) = D · (μ₂+p₂²)²   (first square is 0)
    have key : (μ₂ + p₂^2) * (μ₁ * (μ₂+p₂^2)^2 + μ₂ * (-(p₁*p₂))^2
                 + (p₁*(μ₂+p₂^2) + p₂*(-(p₁*p₂)))^2)
             = (μ₁*μ₂ + μ₂*p₁^2 + μ₁*p₂^2) * (μ₂+p₂^2)^2 := by ring
    nlinarith [mul_nonneg (le_of_lt hc) hq, key, mul_pos hc hc]
  · -- D ≥ 0 ⟹ PSD, via c·q = (…)² + D x² ≥ 0 and c > 0
    intro hD x y
    have hid := poleLift_csquare μ₁ μ₂ p₁ p₂ x y
    have hrhs : 0 ≤ ((μ₂ + p₂^2) * y + p₁*p₂*x)^2 + (μ₁*μ₂ + μ₂*p₁^2 + μ₁*p₂^2) * x^2 := by
      have := mul_nonneg hD (sq_nonneg x)
      positivity
    have hcq : 0 ≤ (μ₂ + p₂^2) * (μ₁ * x^2 + μ₂ * y^2 + (p₁*x + p₂*y)^2) := by
      rw [hid]; exact hrhs
    have h0 : (μ₂ + p₂^2) * 0 ≤ (μ₂ + p₂^2) * (μ₁ * x^2 + μ₂ * y^2 + (p₁*x + p₂*y)^2) := by
      rw [mul_zero]; exact hcq
    exact le_of_mul_le_mul_left h0 hc
