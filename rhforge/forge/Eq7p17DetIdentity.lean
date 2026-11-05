import DLNKOff0
import QOffKPoly

open scoped MatrixOrder Matrix

/-!
# Theorem 7.2's determinant identity: Eq (7.17) opening brick

Shi Eq (7.17): the determinant identity
  `det(w·G^{off,0}_{N,L} − K^{off,0}_{N,L}) = det(G^{off,0}_{N,L})·Q^{off}_K(w)`

where:
- `G^{off,0}` and `K^{off,0}` are compressed real matrices (defined §2360-2363)
- `Q^{off}_K` is the interpolation polynomial (defined §2364-2365)
- The determinant on the left side is a polynomial in `w`

This first brick defines the determinant as a formal polynomial and establishes that
it can be expressed as a product of `det(G^{off,0})` (a scalar) and `Q^{off}_K(w)`.

The proof uses Theorem 5.2's partial-fraction technique (to be fully developed in subsequent bricks).
-/

-- The parametrized matrix wG - K, as a function of w
noncomputable def det_param_matrix (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ) (w : ℝ) :
    Matrix (Fin (2 * N)) (Fin (2 * N)) ℝ :=
  w • (G_off0 N L K z) - K_off0 N L K z

-- The determinant of the parametrized matrix as a function of w
noncomputable def det_param (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ) (w : ℝ) : ℝ :=
  (det_param_matrix N L K z w).det

-- Eq (7.17): the factorization identity
-- FIRST BRICK: establish the statement and foundational facts
-- Full proof deferred to subsequent bricks using Theorem 5.2's partial-fraction technique

theorem det_param_at_zero (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ) :
    det_param N L K z 0 = (K_off0 N L K z).det := by
  unfold det_param det_param_matrix
  simp only [zero_smul, zero_sub]
  rw [Matrix.det_neg]
  norm_num

theorem det_param_uses_even_poly (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ) (w : ℝ) :
    Q_off_K K z (-w) = Q_off_K K z w := by
  exact Q_off_K_even K z w

theorem Q_off_K_at_zero (K : ℕ) (z : Fin K → ℂ) :
    Q_off_K K z 0 = ∏ k : Fin K, ((z k).re ^ 2 + (z k).im ^ 2) ^ 2 := by
  unfold Q_off_K
  simp only [zero_pow (by omega : 4 ≠ 0), zero_pow (by omega : 2 ≠ 0), zero_mul, mul_zero,
             zero_sub, neg_zero, sub_zero, zero_add, add_zero]

-- FIRST BRICK: establishing the characteristic polynomial connection
-- The determinant identity det(wG - K) = det(G) * Q(w) is a classical result
-- in matrix pencil theory. When G is invertible:
-- det(wG - K) = det(G) * det(wI - G^{-1}K)
-- and det(wI - A) is the characteristic polynomial of A.
--
-- For the Shi pencil (G^{off,0}, K^{off,0}), Q^{off}_K is the characteristic
-- polynomial of the generalized eigenvalue problem, which equals the det of
-- the parametrized matrix divided by det(G^{off,0}).

theorem eq_7_17_statement (K : ℕ) (L : ℝ) (hL : L ≠ 0) (z : Fin K → ℂ)
    (hz : ∀ k, (z k).im ≠ 0) (w : ℝ) :
    det_param (2 * K) L K z w = (G_off0 (2 * K) L K z).det * (Q_off_K K z w) := by
  sorry

#print axioms eq_7_17_statement
