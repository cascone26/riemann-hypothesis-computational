import QOffKPoly
import Mathlib

open scoped MatrixOrder Matrix

/-!
# Q_off_K: Real Quadratic Factorization

MATHEMATICAL FOUNDATION (VERIFIED ALGEBRAICALLY):
==================================================

Each quartic factor of Q_off_K can be factored as a product of TWO REAL QUADRATICS:

    w⁴ - 2(γ² - η²)w² + (γ² + η²)² = (w² - 2γw + γ² + η²) · (w² + 2γw + γ² + η²)

where γ = Re(z_k) and η = Im(z_k).

ALGEBRAIC VERIFICATION:
Expanding the right side:
  (w² - 2γw + γ² + η²)(w² + 2γw + γ² + η²)
  = w⁴ + 2γw³ + (γ² + η²)w²
    - 2γw³ - 4γ²w² - 2γ(γ² + η²)w
    + (γ² + η²)w² + 2γ(γ² + η²)w
    + (γ² + η²)²
  = w⁴ + (γ² + η² - 4γ² + γ² + η²)w² + (γ² + η²)²
  = w⁴ + (-2γ² + 2η²)w² + (γ² + η²)²
  = w⁴ - 2(γ² - η²)w² + (γ² + η²)²  ✓

KEY PROPERTIES OF EACH QUADRATIC FACTOR (when η ≠ 0):

For the quadratic q₁(w) = w² - 2γw + γ² + η²:
  - Roots: w = γ ± √(γ² - (γ² + η²)) = γ ± √(-η²) = γ ± iη (complex)
  - Since η ≠ 0, both roots are nonreal
  - Discriminant: (2γ)² - 4(γ² + η²) = 4γ² - 4γ² - 4η² = -4η² < 0
  - Hence q₁(w) > 0 for all real w

For the quadratic q₂(w) = w² + 2γw + γ² + η²:
  - Roots: w = -γ ± √(γ² - (γ² + η²)) = -γ ± iη (complex)
  - Since η ≠ 0, both roots are nonreal
  - Discriminant: (2γ)² - 4(γ² + η²) = -4η² < 0
  - Hence q₂(w) > 0 for all real w

Therefore, each quartic factor is a product of two positive quadratics with no real roots,
confirming the positivity result in QOffKPositivity.lean from a real-factorization perspective.

LEAN FORMALIZATION:
====================
-/

/-- The real quadratic factorization lemma: verify that expanding the product of two real quadratics
    recovers the original quartic. -/
lemma real_quadratic_factorization_identity (w γ η : ℝ) :
    w ^ 4 - 2 * (γ ^ 2 - η ^ 2) * w ^ 2 + (γ ^ 2 + η ^ 2) ^ 2 =
    (w ^ 2 - 2 * γ * w + γ ^ 2 + η ^ 2) * (w ^ 2 + 2 * γ * w + γ ^ 2 + η ^ 2) := by
  ring

/-- Each real quadratic factor has negative discriminant when η ≠ 0, hence nonreal roots. -/
lemma first_quadratic_discriminant (γ η : ℝ) (_hη : η ≠ 0) :
    (2 * γ) ^ 2 - 4 * (γ ^ 2 + η ^ 2) = -4 * η ^ 2 := by
  ring

/-- The first quadratic w² - 2γw + γ² + η² is positive for all real w when η ≠ 0. -/
lemma first_quadratic_pos (w γ η : ℝ) (hη : η ≠ 0) :
    0 < w ^ 2 - 2 * γ * w + γ ^ 2 + η ^ 2 := by
  have h : w ^ 2 - 2 * γ * w + γ ^ 2 + η ^ 2 = (w - γ) ^ 2 + η ^ 2 := by ring
  rw [h]
  have h_sq_nonneg : 0 ≤ (w - γ) ^ 2 := sq_nonneg _
  have h_eta_sq_pos : 0 < η ^ 2 := sq_pos_of_ne_zero hη
  linarith

/-- The second quadratic w² + 2γw + γ² + η² is positive for all real w when η ≠ 0. -/
lemma second_quadratic_pos (w γ η : ℝ) (hη : η ≠ 0) :
    0 < w ^ 2 + 2 * γ * w + γ ^ 2 + η ^ 2 := by
  have h : w ^ 2 + 2 * γ * w + γ ^ 2 + η ^ 2 = (w + γ) ^ 2 + η ^ 2 := by ring
  rw [h]
  have h_sq_nonneg : 0 ≤ (w + γ) ^ 2 := sq_nonneg _
  have h_eta_sq_pos : 0 < η ^ 2 := sq_pos_of_ne_zero hη
  linarith

/-- Both quadratic factors are nonzero for all real w when η ≠ 0. -/
lemma first_quadratic_ne_zero (w γ η : ℝ) (hη : η ≠ 0) :
    w ^ 2 - 2 * γ * w + γ ^ 2 + η ^ 2 ≠ 0 :=
  (first_quadratic_pos w γ η hη).ne'

lemma second_quadratic_ne_zero (w γ η : ℝ) (hη : η ≠ 0) :
    w ^ 2 + 2 * γ * w + γ ^ 2 + η ^ 2 ≠ 0 :=
  (second_quadratic_pos w γ η hη).ne'

/-- Product of the two quadratics is positive when η ≠ 0. -/
lemma quadratic_product_pos (w γ η : ℝ) (hη : η ≠ 0) :
    0 < (w ^ 2 - 2 * γ * w + γ ^ 2 + η ^ 2) * (w ^ 2 + 2 * γ * w + γ ^ 2 + η ^ 2) :=
  mul_pos (first_quadratic_pos w γ η hη) (second_quadratic_pos w γ η hη)

/-- The quartic factor, expressed as a product of two real quadratics, is positive when η ≠ 0. -/
theorem quartic_as_real_quadratic_product (w γ η : ℝ) (hη : η ≠ 0) :
    0 < w ^ 4 - 2 * (γ ^ 2 - η ^ 2) * w ^ 2 + (γ ^ 2 + η ^ 2) ^ 2 := by
  rw [real_quadratic_factorization_identity]
  exact quadratic_product_pos w γ η hη

/-- Formalized statement: Q_off_K can be viewed as a product of real quadratic factors. -/
theorem Q_off_K_real_quadratic_factorization (K : ℕ) (z : Fin K → ℂ) (w : ℝ) :
    Q_off_K K z w = ∏ k : Fin K,
      ((w ^ 2 - 2 * (z k).re * w + (z k).re ^ 2 + (z k).im ^ 2) *
       (w ^ 2 + 2 * (z k).re * w + (z k).re ^ 2 + (z k).im ^ 2)) := by
  unfold Q_off_K
  apply Finset.prod_congr rfl
  intro k _
  exact real_quadratic_factorization_identity w (z k).re (z k).im

#print axioms real_quadratic_factorization_identity
#print axioms first_quadratic_pos
#print axioms second_quadratic_pos
#print axioms quartic_as_real_quadratic_product
#print axioms Q_off_K_real_quadratic_factorization
