import Mathlib

open Matrix

namespace CMKeyLemma

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **Matrix determinant lemma (rank-one, outer-product form).**

For an invertible real matrix `A` and vectors `u, v`,
    `det (A + |u⟩⟨v|) = det A · (1 + ⟨v, A⁻¹ u⟩)`,
with `|u⟩⟨v| = vecMulVec u v`.  This is the algebraic engine of CM 2025 Lemma "key" part (iii):
applied with `A = D − s` and `u = −Dξ`, `v = η`, it yields
`Det(D'−s) = Det(D−s)·(1 − ⟨η,(D−s)⁻¹Dξ⟩)`, the spectral-realization determinant.  Machine-verified,
axiom-clean. -/
theorem det_add_vecMulVec (A : Matrix n n ℝ) (u v : n → ℝ) (hA : IsUnit A.det) :
    (A + Matrix.vecMulVec u v).det = A.det * (1 + v ⬝ᵥ (A⁻¹ *ᵥ u)) := by
  have hfac : A * (1 + A⁻¹ * Matrix.vecMulVec u v) = A + Matrix.vecMulVec u v := by
    rw [Matrix.mul_add, Matrix.mul_one, ← Matrix.mul_assoc, Matrix.mul_nonsing_inv A hA,
      Matrix.one_mul]
  rw [← hfac, Matrix.det_mul, Matrix.mul_vecMulVec, Matrix.vecMulVec_eq (ι := Unit),
    Matrix.det_one_add_replicateCol_mul_replicateRow]

#print axioms det_add_vecMulVec

end CMKeyLemma
