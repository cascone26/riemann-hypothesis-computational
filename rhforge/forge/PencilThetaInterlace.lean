/-
PencilThetaInterlace.lean (HP RH center, LEDGER S2903, 2026-11-07)
S2884 composition step: if B_{N+1} = [[B_N, 0],[yᵀ, β]] (top-right column zero — the lower-triangularity of Shi's
unshifted contrast pencil), then the eigenvalues θ² of B_N B_Nᵀ interlace with those of B_{N+1} B_{N+1}ᵀ:
  θ²_{k+1}(N+1) ≤ θ²_k(N) ≤ θ²_k(N+1).
Composes BlockLowerGram (S2900) with MatrixCauchyInterlacing (S2902).
-/
import Mathlib
import BlockLowerGram
import MatrixCauchyInterlacing

open Matrix Module

namespace ThetaInterlace

/-- eigenvalues depend only on the matrix (proof irrelevance + subst) -/
lemma eigenvalues_congr {n : ℕ} {X Y : Matrix (Fin n) (Fin n) ℝ} (h : X = Y)
    (hX : (toEuclideanLin X).IsSymmetric) (hY : (toEuclideanLin Y).IsSymmetric)
    (hn : finrank ℝ (EuclideanSpace ℝ (Fin n)) = n) :
    hX.eigenvalues hn = hY.eigenvalues hn := by
  subst h; rfl

theorem theta_sq_interlace {m : ℕ} (B : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ)
    (hcol : ∀ i : Fin m, B i.castSucc (Fin.last m) = 0)
    (hBig : (toEuclideanLin (B * Bᵀ)).IsSymmetric)
    (hSmall : (toEuclideanLin (B.submatrix Fin.castSucc Fin.castSucc *
      (B.submatrix Fin.castSucc Fin.castSucc)ᵀ)).IsSymmetric) (k : Fin m) :
    hBig.eigenvalues (MatrixInterlace.finrank_e (m + 1)) k.succ
        ≤ hSmall.eigenvalues (MatrixInterlace.finrank_e m) k ∧
      hSmall.eigenvalues (MatrixInterlace.finrank_e m) k
        ≤ hBig.eigenvalues (MatrixInterlace.finrank_e (m + 1)) k.castSucc := by
  have hA : (B * Bᵀ).IsHermitian := by
    have := isHermitian_mul_conjTranspose_self B
    rwa [conjTranspose_eq_transpose_of_trivial] at this
  obtain ⟨h1, h2⟩ := MatrixInterlace.matrix_interlace (B * Bᵀ) hA k
  have e := BlockGram.gram_principal_block B hcol
  have hsub : (toEuclideanLin ((B * Bᵀ).submatrix Fin.castSucc Fin.castSucc)).IsSymmetric :=
    (Matrix.isHermitian_iff_isSymmetric.mp (hA.submatrix Fin.castSucc))
  have eq := eigenvalues_congr e hsub hSmall (MatrixInterlace.finrank_e m)
  rw [eq] at h1 h2
  exact ⟨h1, h2⟩

end ThetaInterlace

#print axioms ThetaInterlace.theta_sq_interlace
