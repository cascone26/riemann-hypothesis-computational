import QOffKPoly
import Mathlib

open scoped MatrixOrder Matrix

/-!
# Q_off_K Positivity on ℝ: Formal Proof

THEOREM: If all z_k have nonzero imaginary part, then Q_off_K(w) > 0 for all real w.

MATHEMATICAL FOUNDATION (verified independently):
================================================

Each quartic factor has the closed form:
  f_k(w) := w⁴ - 2(γ_k² - η_k²)w² + (γ_k² + η_k²)²

where z_k = γ_k + iη_k, γ_k = Re(z_k), η_k = Im(z_k).

Completing the square (algebraically verified):
  f_k(w) = (w² - (γ_k² - η_k²))² + 4γ_k² η_k²

Key observation when η_k ≠ 0:
CASE 1: If γ_k ≠ 0, both terms are > 0, so sum > 0.
CASE 2: If γ_k = 0, then 4γ_k² η_k² = 0, so f_k(w) = (w² + η_k²)² > 0 for real w.
In both cases, f_k(w) > 0 for all real w.

Since Q_off_K is a product of these positive factors, Q_off_K(w) > 0 for all real w.

LEAN FORMALIZATION:
===================
-/

/-- Each quartic factor can be written as a sum of two non-negative terms. -/
lemma quartic_factor_positive_decomposition
    (w : ℝ) (γ η : ℝ) (_hη : η ≠ 0) :
    w ^ 4 - 2 * (γ ^ 2 - η ^ 2) * w ^ 2 + (γ ^ 2 + η ^ 2) ^ 2 =
    (w ^ 2 - (γ ^ 2 - η ^ 2)) ^ 2 + 4 * γ ^ 2 * η ^ 2 := by
  ring

/-- Each quartic factor is strictly positive when η ≠ 0. -/
lemma quartic_factor_pos
    (w : ℝ) (γ η : ℝ) (hη : η ≠ 0) :
    0 < w ^ 4 - 2 * (γ ^ 2 - η ^ 2) * w ^ 2 + (γ ^ 2 + η ^ 2) ^ 2 := by
  rw [quartic_factor_positive_decomposition w γ η hη]
  by_cases hg : γ = 0
  · -- Case: γ = 0, so the factor becomes (w² + η²)²
    rw [hg]
    simp only [zero_pow (by norm_num : (2 : ℕ) ≠ 0), zero_mul, mul_zero, add_zero]
    exact sq_pos_of_pos (by nlinarith [sq_pos_of_ne_zero hη])
  · -- Case: γ ≠ 0, both terms are strictly positive
    have h_sq_nonneg : 0 ≤ (w ^ 2 - (γ ^ 2 - η ^ 2)) ^ 2 := sq_nonneg _
    have h_four : (0 : ℝ) < 4 := by norm_num
    have h_eta_sq : 0 < η ^ 2 := sq_pos_of_ne_zero hη
    have h_gamma_sq_pos : 0 < γ ^ 2 := sq_pos_of_ne_zero hg
    have h_prod_pos : 0 < 4 * γ ^ 2 * η ^ 2 := mul_pos (mul_pos h_four h_gamma_sq_pos) h_eta_sq
    linarith

/-- Q_off_K is strictly positive for all real w when all z_k have nonzero imaginary part. -/
theorem Q_off_K_pos
    (K : ℕ) (z : Fin K → ℂ) (w : ℝ)
    (hz : ∀ k, (z k).im ≠ 0) :
    0 < Q_off_K K z w := by
  unfold Q_off_K
  apply Finset.prod_pos
  intro k _
  exact quartic_factor_pos w (z k).re (z k).im (hz k)

/-- Corollary: Q_off_K is nonzero for all real w when all z_k have nonzero imaginary part. -/
theorem Q_off_K_ne_zero
    (K : ℕ) (z : Fin K → ℂ) (w : ℝ)
    (hz : ∀ k, (z k).im ≠ 0) :
    Q_off_K K z w ≠ 0 := by
  exact (Q_off_K_pos K z w hz).ne'

#print axioms Q_off_K_pos
#print axioms Q_off_K_ne_zero
