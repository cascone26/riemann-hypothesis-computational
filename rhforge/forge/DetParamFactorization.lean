import Eq7p17DetIdentity
import Mathlib

open Matrix

/-!
# Eq (7.17) — Determinant Factorization Approach

Direct computation strategy: factorize `det(w·G_off0 - K_off0)` using the standard
determinant identity for invertible matrices:

If G is invertible: `det(w·G - K) = det(G) · det(w·I - G⁻¹·K)`

This reduces Eq (7.17) to a characteristic polynomial problem.

KEY FACT: If G is invertible, then det(w·G - K) factors as det(G) times the characteristic
polynomial of G⁻¹·K.
-/

-- Standard determinant factorization for invertible matrices
-- For invertible G: det(w·G - K) = det(G) · det(w·I - G⁻¹·K)
-- This follows from the matrix identity w·G - K = G·(w·I - G⁻¹·K)
-- and the multiplicativity of the determinant
theorem matrix_wG_sub_K_eq_mul {n : Type*} [Fintype n] [DecidableEq n]
    (G K : Matrix n n ℝ) (w : ℝ) (hG : IsUnit G.det) :
    w • G - K = G * (w • 1 - G⁻¹ * K) := by
  rw [Matrix.mul_sub, Matrix.mul_smul, Matrix.mul_one, mul_nonsing_inv_cancel_left G K hG]

theorem det_wG_sub_K_factorization {n : Type*} [Fintype n] [DecidableEq n]
    (G K : Matrix n n ℝ) (w : ℝ) (hG : IsUnit G.det) :
    (w • G - K).det = G.det * (w • 1 - G⁻¹ * K).det := by
  rw [matrix_wG_sub_K_eq_mul G K w hG, Matrix.det_mul]

-- Apply to our concrete case: det(w·G_off0 - K_off0) = det(G_off0) * det(w·1 - G_off0⁻¹·K_off0)
theorem det_param_factorization (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ) (w : ℝ)
    (hG : IsUnit (G_off0 N L K z).det) :
    det_param N L K z w = (G_off0 N L K z).det * ((w • 1 - (G_off0 N L K z)⁻¹ * (K_off0 N L K z)).det) := by
  unfold det_param det_param_matrix
  exact det_wG_sub_K_factorization (G_off0 N L K z) (K_off0 N L K z) w hG

-- The reduction: Shows that Eq (7.17) follows IF:
-- (1) G_off0 is invertible
-- (2) The characteristic polynomial of G_off0⁻¹·K_off0 equals Q_off_K(w)
theorem eq_7_17_reduction (N : ℕ) (L : ℝ) (_hL : L ≠ 0) (K : ℕ) (z : Fin K → ℂ)
    (_hz : ∀ k, (z k).im ≠ 0) (w : ℝ)
    (hG : IsUnit (G_off0 N L K z).det)
    (h_char : (w • 1 - (G_off0 N L K z)⁻¹ * (K_off0 N L K z)).det = Q_off_K K z w) :
    det_param N L K z w = (G_off0 N L K z).det * (Q_off_K K z w) := by
  rw [← h_char]
  exact det_param_factorization N L K z w hG

#print axioms det_wG_sub_K_factorization
#print axioms det_param_factorization
#print axioms eq_7_17_reduction
