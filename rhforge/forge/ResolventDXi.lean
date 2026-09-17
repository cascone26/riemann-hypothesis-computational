import Mathlib

open Matrix

namespace CMKeyLemma

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **Resolvent identity for CM Lemma "key" part (iii).**

With `A = D − s·1` invertible, `Dξ = (D−s)ξ + sξ` gives
    `(D − s·1)⁻¹ · (D · ξ) = ξ + s · (D − s·1)⁻¹ · ξ`.
This is the second half of CM's determinant computation (iii): fed into the matrix determinant
lemma `det_add_vecMulVec`, it turns `1 − ⟨η,(D−s)⁻¹Dξ⟩` into `−s⟨η,(D−s)⁻¹ξ⟩` (using `⟨η,ξ⟩=1`).
Machine-verified, axiom-clean. -/
theorem resolvent_D_mulVec (D : Matrix n n ℝ) (ξ : n → ℝ) (s : ℝ)
    (hA : IsUnit (D - s • (1 : Matrix n n ℝ)).det) :
    (D - s • (1 : Matrix n n ℝ))⁻¹ *ᵥ (D *ᵥ ξ)
      = ξ + s • ((D - s • (1 : Matrix n n ℝ))⁻¹ *ᵥ ξ) := by
  set A := D - s • (1 : Matrix n n ℝ) with hAdef
  have hDdecomp : D *ᵥ ξ = A *ᵥ ξ + s • ξ := by
    rw [hAdef, Matrix.sub_mulVec, Matrix.smul_mulVec, Matrix.one_mulVec]; abel
  rw [hDdecomp, Matrix.mulVec_add, Matrix.mulVec_smul]
  congr 1
  rw [Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul A hA, Matrix.one_mulVec]

#print axioms resolvent_D_mulVec

end CMKeyLemma
