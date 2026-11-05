import GOff0Symm
import DLNKOff0

open Matrix

/-!
# Off-critical-line lane: `G_off0` Invertibility Investigation

Investigation of whether `G_off0 N L K z` is always invertible for generic z with nonzero imaginary parts.

G_off0 is defined as: G_off0 = C^T * S_T_off_real * C
where:
- C = C_matrix (2 * N + 1) : Matrix (Fin(2N+1)) (Fin(2N))
- S_T_off_real : Matrix (Fin(2N+1)) (Fin(2N+1)) is a real symmetric matrix, generally indefinite

The structure is a congruence of an indefinite symmetric matrix.

Key facts:
1. S_T_off_real is symmetric (proven)
2. S_T_off_real has a null vector e_off (from STOffFinalAssembly.lean)
3. C_matrix is injective in its mulVec action (from ContrastMatrixPacked comments)

Hypothesis: G_off0 is invertible because the null vector e_off of S_T_off_real
is orthogonal to the image of C_matrix (i.e., not in the image of C).
-/

-- Attempted theorem: G_off0 is invertible
-- This would require showing that the kernel of G_off0 is trivial
theorem G_off0_isUnit (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ) (hL : L ≠ 0) (hz : ∀ k, (z k).im ≠ 0) :
    IsUnit (G_off0 N L K z).det := by
  -- Strategy: show that G_off0 has full rank
  -- For a congruence G = C^T * S * C, if C has full column rank and S is nonsingular on
  -- the complement of its null space, then G is invertible iff the null space of S
  -- is orthogonal to the image of C.
  sorry

#print axioms G_off0_isUnit
