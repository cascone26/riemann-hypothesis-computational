/-
BlockLowerGram.lean (HP RH center, LEDGER S2900, 2026-11-07)
S2884 structural step: if B_{N+1} = [[B_N, 0],[yᵀ, β]] (top-right column zero — the lower-triangularity of
Shi's unshifted contrast pencil, proved on paper in S2884, core algebra in Lean S2885), then B_N B_Nᵀ is the
leading principal block of B_{N+1} B_{N+1}ᵀ. With Cauchy interlacing (S2899) this gives one-step interlacing
of the pencil's singular values θ (θ² = eigenvalues of B Bᵀ).
-/
import Mathlib

open Matrix

namespace BlockGram

theorem gram_principal_block {m : ℕ} (B : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ)
    (hcol : ∀ i : Fin m, B i.castSucc (Fin.last m) = 0) :
    (B * Bᵀ).submatrix Fin.castSucc Fin.castSucc
      = B.submatrix Fin.castSucc Fin.castSucc * (B.submatrix Fin.castSucc Fin.castSucc)ᵀ := by
  ext i j
  simp only [submatrix_apply, mul_apply, transpose_apply]
  rw [Fin.sum_univ_castSucc]
  simp [hcol i]

/-- Same statement with an explicit block decomposition via `Matrix.fromBlocks`-free indexing:
the Gram matrix of the restriction equals the restriction of the Gram matrix. -/
theorem gram_principal_block_symm {m : ℕ} (B : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ)
    (hcol : ∀ i : Fin m, B i.castSucc (Fin.last m) = 0) :
    ((B * Bᵀ).submatrix Fin.castSucc Fin.castSucc).IsSymm := by
  rw [gram_principal_block B hcol]
  exact isSymm_mul_transpose_self _

end BlockGram

#print axioms BlockGram.gram_principal_block
#print axioms BlockGram.gram_principal_block_symm
