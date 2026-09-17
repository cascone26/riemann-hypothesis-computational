import DetAddVecMulVec
import ResolventDXi

open Matrix

namespace CMKeyLemma

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **CM 2025 Lemma "key" part (iii), core determinant identity.**

For the perturbed operator `D' = D − |Dξ⟩⟨η|` (`= D − vecMulVec (Dξ) η`), with `A = D − s·1`
invertible and the normalization `⟨η,ξ⟩ = 1`,

    `det(D' − s·1) = −s · det(D − s·1) · ⟨η, (D − s·1)⁻¹ ξ⟩`.

Composed from the matrix determinant lemma (`det_add_vecMulVec`) and the resolvent identity
(`resolvent_D_mulVec`).  This is the spectral-realization determinant: the eigenvalues of `D''`
(the zeros of the RHS at `s ∉ spec D`) are the zeros of `s ↦ ⟨η,(D−s)⁻¹ξ⟩ = Σ_j (j−s)⁻¹ξ_j = ξ̂(s)`.
Machine-verified, axiom-clean. -/
theorem detDprime_core (D : Matrix n n ℝ) (η ξ : n → ℝ) (s : ℝ)
    (hA : IsUnit (D - s • (1 : Matrix n n ℝ)).det) (hη : η ⬝ᵥ ξ = 1) :
    (D - Matrix.vecMulVec (D *ᵥ ξ) η - s • (1 : Matrix n n ℝ)).det
      = -s * (D - s • (1 : Matrix n n ℝ)).det
          * (η ⬝ᵥ ((D - s • (1 : Matrix n n ℝ))⁻¹ *ᵥ ξ)) := by
  set A := D - s • (1 : Matrix n n ℝ) with hAdef
  have hres := resolvent_D_mulVec D ξ s hA
  rw [← hAdef] at hres
  have hstep : D - Matrix.vecMulVec (D *ᵥ ξ) η - s • (1 : Matrix n n ℝ)
      = A + Matrix.vecMulVec (-(D *ᵥ ξ)) η := by
    simp only [hAdef, Matrix.neg_vecMulVec]; abel
  rw [hstep, det_add_vecMulVec A (-(D *ᵥ ξ)) η hA, Matrix.mulVec_neg, dotProduct_neg, hres,
    dotProduct_add, dotProduct_smul, hη]
  simp only [smul_eq_mul]; ring

#print axioms detDprime_core

end CMKeyLemma
