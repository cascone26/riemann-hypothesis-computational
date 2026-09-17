import Mathlib

open Matrix

namespace CMKeyLemma

variable {n : Type*} [Fintype n] [DecidableEq n]

private lemma vmv_act (a b w : n → ℝ) :
    Matrix.vecMulVec a b *ᵥ w = (b ⬝ᵥ w) • a := by
  funext i
  simp only [Matrix.mulVec, Matrix.vecMulVec_apply, dotProduct, Pi.smul_apply, smul_eq_mul,
    Finset.sum_mul]
  exact Finset.sum_congr rfl (fun j _ => by ring)

/-- **CM 2025 Lemma "key": the kernel fact `D'ξ = 0`.**

The perturbed operator `D' = D − |Dξ⟩⟨η|` annihilates `ξ` whenever `⟨η,ξ⟩ = 1`:
    `(D − vecMulVec (Dξ) η) *ᵥ ξ = 0`.
CM use exactly this ("by construction `D'ξ=0`, so `D'` induces `D''` on the quotient `E_N/ℂξ`"):
it is the foundation both of the self-adjoint quotient `D''` (part ii) and of the determinant
factorization `Det(D'−s) = −s·Det(D''−s)` (the `det0` step of part iii).  Unconditional — needs only
the normalization, not even-simplicity.  Machine-verified, axiom-clean. -/
theorem Dprime_mulVec_xi (D : Matrix n n ℝ) (η ξ : n → ℝ) (hη : η ⬝ᵥ ξ = 1) :
    (D - Matrix.vecMulVec (D *ᵥ ξ) η) *ᵥ ξ = 0 := by
  rw [Matrix.sub_mulVec, vmv_act, hη, one_smul, sub_self]

#print axioms Dprime_mulVec_xi

end CMKeyLemma
