/-
MatrixCauchyInterlacing.lean (HP RH center, LEDGER S2902, 2026-11-07)
Matrix form of Cauchy interlacing: a real symmetric (m+1)×(m+1) matrix A and its leading m×m principal
submatrix A' = A.submatrix castSucc castSucc have interlacing eigenvalues λ_{k+1}(A) ≤ μ_k(A') ≤ λ_k(A)
(decreasing order, as operators on EuclideanSpace). Via CauchyInterlaceIso (S2901) with ι = pad-by-zero.
-/
import Mathlib
import CauchyInterlacingIso

open Matrix Module

namespace MatrixInterlace

variable {m : ℕ}

/-- pad a vector of length m with a zero at the last index -/
def padFun (v : Fin m → ℝ) : Fin (m + 1) → ℝ := Fin.lastCases 0 v

lemma padFun_castSucc (v : Fin m → ℝ) (j : Fin m) : padFun v j.castSucc = v j := by
  simp [padFun]
lemma padFun_last (v : Fin m → ℝ) : padFun v (Fin.last m) = 0 := by simp [padFun]

noncomputable def padLin : EuclideanSpace ℝ (Fin m) →ₗ[ℝ] EuclideanSpace ℝ (Fin (m + 1)) where
  toFun v := (WithLp.equiv 2 _).symm (padFun ((WithLp.equiv 2 _) v))
  map_add' v w := by
    ext i; refine Fin.lastCases ?_ (fun j => ?_) i <;> simp [padFun]
  map_smul' c v := by
    ext i; refine Fin.lastCases ?_ (fun j => ?_) i <;> simp [padFun]

lemma padLin_apply_castSucc (v : EuclideanSpace ℝ (Fin m)) (j : Fin m) :
    padLin v j.castSucc = v j := by simp [padLin, padFun]
lemma padLin_apply_last (v : EuclideanSpace ℝ (Fin m)) : padLin v (Fin.last m) = 0 := by
  simp [padLin, padFun]

noncomputable def padIso : EuclideanSpace ℝ (Fin m) →ₗᵢ[ℝ] EuclideanSpace ℝ (Fin (m + 1)) where
  toLinearMap := padLin
  norm_map' v := by
    rw [EuclideanSpace.norm_eq, EuclideanSpace.norm_eq, Fin.sum_univ_castSucc]
    simp [padLin_apply_castSucc, padLin_apply_last]

lemma inner_eq_sum (x y : EuclideanSpace ℝ (Fin m)) : @inner ℝ _ _ x y = ∑ i, x i * y i := by
  simp [PiLp.inner_apply, mul_comm]

theorem compression (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (v w : EuclideanSpace ℝ (Fin m)) :
    @inner ℝ _ _ (toEuclideanLin (A.submatrix Fin.castSucc Fin.castSucc) v) w
      = @inner ℝ _ _ (toEuclideanLin A (padIso v)) (padIso w) := by
  rw [toEuclideanLin_apply, toEuclideanLin_apply]
  simp only [PiLp.inner_apply]
  rw [Fin.sum_univ_castSucc (fun i => _)]
  simp [padIso, padLin_apply_castSucc, padLin_apply_last, mulVec, dotProduct,
    Fin.sum_univ_castSucc, submatrix_apply]

lemma finrank_e (n : ℕ) : finrank ℝ (EuclideanSpace ℝ (Fin n)) = n := finrank_euclideanSpace_fin

theorem matrix_interlace (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (hA : A.IsHermitian) (k : Fin m) :
    let hT := (isHermitian_iff_isSymmetric.mp hA)
    let hS := (isHermitian_iff_isSymmetric.mp (hA.submatrix Fin.castSucc))
    hT.eigenvalues (finrank_e (m + 1)) k.succ ≤ hS.eigenvalues (finrank_e m) k ∧
      hS.eigenvalues (finrank_e m) k ≤ hT.eigenvalues (finrank_e (m + 1)) k.castSucc := by
  intro hT hS
  have hST := compression A
  exact ⟨CauchyInterlaceIso.interlace_lower padIso hT hS hST (finrank_e (m + 1)) (finrank_e m) k,
    CauchyInterlaceIso.interlace_upper padIso hT hS hST (finrank_e (m + 1)) (finrank_e m) k⟩

end MatrixInterlace

#print axioms MatrixInterlace.matrix_interlace
