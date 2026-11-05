import Mathlib

open Matrix

lemma matrix_identity_for_factorization {n : Type*} [Fintype n] [DecidableEq n]
    (G K : Matrix n n ℝ) (w : ℝ) (hG : IsUnit G.det) :
    w • G - K = G * (w • 1 - G⁻¹ * K) := by
  rw [Matrix.mul_sub, Matrix.mul_smul, Matrix.mul_one, mul_nonsing_inv_cancel_left G K hG]

theorem det_wG_sub_K_factorization {n : Type*} [Fintype n] [DecidableEq n]
    (G K : Matrix n n ℝ) (w : ℝ) (hG : IsUnit G.det) :
    (w • G - K).det = G.det * (w • 1 - G⁻¹ * K).det := by
  rw [matrix_identity_for_factorization G K w hG, Matrix.det_mul]

#print axioms matrix_identity_for_factorization
#print axioms det_wG_sub_K_factorization
