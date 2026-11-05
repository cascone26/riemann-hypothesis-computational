/-
PencilDisplacementSwap.lean (HP RH center, LEDGER S2885, 2026-11-05)
Core algebra of S2884's triangularity lemma for Shi's unshifted contrast pencil:
if DS − SD = b1ᵀ − 1bᵀ (Shi Eq 3.31, rank-2 displacement) then on zero-mean vectors the
differentiation form uᵀ(SD)v equals (Dᵀu)ᵀSv; hence it vanishes whenever Dᵀu lies in a subspace
that v is S-orthogonal to.
-/
import Mathlib

open Matrix

namespace PencilSwap

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- displacement swap on zero-mean vectors -/
theorem displacement_swap (S D : Matrix n n ℝ) (b : n → ℝ)
    (hdisp : D * S - S * D = vecMulVec b (fun _ => 1) - vecMulVec (fun _ => 1) b)
    (u v : n → ℝ) (hu : ∑ i, u i = 0) (hv : ∑ i, v i = 0) :
    u ⬝ᵥ ((S * D) *ᵥ v) = (Dᵀ *ᵥ u) ⬝ᵥ (S *ᵥ v) := by
  have h1 : (Dᵀ *ᵥ u) ⬝ᵥ (S *ᵥ v) = u ⬝ᵥ ((D * S) *ᵥ v) := by
    rw [← Matrix.vecMul_transpose, transpose_transpose, ← Matrix.dotProduct_mulVec, Matrix.mulVec_mulVec]
  have hS : S * D = D * S - (vecMulVec b (fun _ => 1) - vecMulVec (fun _ => 1) b) := by
    rw [← hdisp]; abel
  have h1v : (fun _ => (1:ℝ)) ⬝ᵥ v = 0 := by simpa [dotProduct] using hv
  have h1u : u ⬝ᵥ (fun _ => (1:ℝ)) = 0 := by simpa [dotProduct] using hu
  rw [h1, hS, Matrix.sub_mulVec, Matrix.sub_mulVec, Matrix.vecMulVec_mulVec, Matrix.vecMulVec_mulVec,
    dotProduct_sub, dotProduct_sub]
  simp only [dotProduct_smul, smul_eq_mul]
  rw [h1v]
  have : u ⬝ᵥ (((b ⬝ᵥ v)) • (fun _ => (1:ℝ))) = (b ⬝ᵥ v) * (u ⬝ᵥ (fun _ => (1:ℝ))) := by
    rw [dotProduct_smul, smul_eq_mul]
  simp [this, h1u]

/-- vanishing corollary: if Dᵀu ∈ W and v is S-orthogonal to W, the differentiation form vanishes -/
theorem form_vanishes (S D : Matrix n n ℝ) (b : n → ℝ)
    (hdisp : D * S - S * D = vecMulVec b (fun _ => 1) - vecMulVec (fun _ => 1) b)
    (W : Submodule ℝ (n → ℝ)) (u v : n → ℝ) (hu : ∑ i, u i = 0) (hv : ∑ i, v i = 0)
    (hDu : Dᵀ *ᵥ u ∈ W) (hvW : ∀ w ∈ W, w ⬝ᵥ (S *ᵥ v) = 0) :
    u ⬝ᵥ ((S * D) *ᵥ v) = 0 := by
  rw [displacement_swap S D b hdisp u v hu hv]
  exact hvW _ hDu

end PencilSwap

#print axioms PencilSwap.displacement_swap
#print axioms PencilSwap.form_vanishes
