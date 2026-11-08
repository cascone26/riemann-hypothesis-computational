/-
MatrixInterlaceAt.lean (HP RH center, LEDGER S2910, 2026-11-08)
Matrix Cauchy interlacing for deleting ANY index p (A' = A.submatrix (succAbove p) (succAbove p)), via the isometry
"insert 0 at position p". Gives λ_min(A) ≤ λ_min(A') for every one-index deletion — so the symmetric growth of the
Weil matrices τ_{λ,N} ⊂ τ_{λ,N+1} (indices −N..N inside −N−1..N+1: delete first, then last) yields ε_{λ,N+1} ≤ ε_{λ,N}.
-/
import Mathlib
import CauchyInterlacingIso
import MatrixCauchyInterlacing

open Matrix Module

namespace InterlaceAt

variable {m : ℕ}

noncomputable def padAtLin (p : Fin (m + 1)) : EuclideanSpace ℝ (Fin m) →ₗ[ℝ] EuclideanSpace ℝ (Fin (m + 1)) where
  toFun v := (WithLp.equiv 2 _).symm (Fin.insertNth p 0 ((WithLp.equiv 2 _) v))
  map_add' v w := by
    ext i; refine Fin.succAboveCases p ?_ (fun j => ?_) i <;> simp
  map_smul' c v := by
    ext i; refine Fin.succAboveCases p ?_ (fun j => ?_) i <;> simp

lemma padAt_succAbove (p : Fin (m + 1)) (v : EuclideanSpace ℝ (Fin m)) (j : Fin m) :
    padAtLin p v (p.succAbove j) = v j := by simp [padAtLin]
lemma padAt_p (p : Fin (m + 1)) (v : EuclideanSpace ℝ (Fin m)) : padAtLin p v p = 0 := by
  simp [padAtLin]

noncomputable def padAtIso (p : Fin (m + 1)) : EuclideanSpace ℝ (Fin m) →ₗᵢ[ℝ] EuclideanSpace ℝ (Fin (m + 1)) where
  toLinearMap := padAtLin p
  norm_map' v := by
    rw [EuclideanSpace.norm_eq, EuclideanSpace.norm_eq, Fin.sum_univ_succAbove _ p]
    simp [padAt_succAbove, padAt_p]

theorem compressionAt (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (p : Fin (m + 1)) (v w : EuclideanSpace ℝ (Fin m)) :
    @inner ℝ _ _ (toEuclideanLin (A.submatrix p.succAbove p.succAbove) v) w
      = @inner ℝ _ _ (toEuclideanLin A (padAtIso p v)) (padAtIso p w) := by
  rw [toEuclideanLin_apply, toEuclideanLin_apply]
  simp only [PiLp.inner_apply]
  rw [Fin.sum_univ_succAbove (fun i => _) p]
  simp [padAtIso, padAt_succAbove, padAt_p, mulVec, dotProduct, Fin.sum_univ_succAbove _ p, submatrix_apply]

theorem min_eig_le_delete (A : Matrix (Fin (m + 2)) (Fin (m + 2)) ℝ) (hA : A.IsHermitian) (p : Fin (m + 2)) :
    (isHermitian_iff_isSymmetric.mp hA).eigenvalues (MatrixInterlace.finrank_e (m + 2)) (Fin.last (m + 1))
      ≤ (isHermitian_iff_isSymmetric.mp (hA.submatrix p.succAbove)).eigenvalues
          (MatrixInterlace.finrank_e (m + 1)) (Fin.last m) := by
  have h := CauchyInterlaceIso.interlace_lower (padAtIso p) (isHermitian_iff_isSymmetric.mp hA)
    (isHermitian_iff_isSymmetric.mp (hA.submatrix p.succAbove)) (compressionAt A p)
    (MatrixInterlace.finrank_e (m + 2)) (MatrixInterlace.finrank_e (m + 1)) (Fin.last m)
  simpa [Fin.succ_last] using h

end InterlaceAt

#print axioms InterlaceAt.min_eig_le_delete
