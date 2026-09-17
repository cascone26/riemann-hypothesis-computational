import Mathlib

open Matrix

namespace CMKeyLemma

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- Action of an outer product on a vector (over ℝ, commutative — avoids the `MulOpposite`
in Mathlib's `vecMulVec_mulVec`). -/
private lemma vecMulVec_act (a b w : n → ℝ) :
    Matrix.vecMulVec a b *ᵥ w = (b ⬝ᵥ w) • a := by
  funext i
  simp only [Matrix.mulVec, Matrix.vecMulVec_apply, dotProduct, Pi.smul_apply, smul_eq_mul,
    Finset.sum_mul]
  exact Finset.sum_congr rfl (fun j _ => by ring)

/-- **CM 2025 (ccm_2511.22755) Lemma "key" part (i), abstract form.**

Setup of CM §835-870: on a finite-dim real space, `D`,`T` are operators whose commutator is the
rank-two form `D·T − T·D = |β⟩⟨η| − |η⟩⟨β|` (their eq. (QD)).  For the minimal (kernel) vector `ξ`
of the even-simple positive form `T` — so `T·ξ = 0` — normalized by `⟨η,ξ⟩ = 1`, and with `⟨β,ξ⟩ = 0`
(β is odd, ξ even, the γ-eigenspaces are orthogonal), one has

    `T · (D · ξ) = −β`.

This is the first step of CM's construction of the perturbed self-adjoint operator `D''` realizing
the ζ-spectrum. Pure finite-dimensional linear algebra; machine-verified, axiom-clean. -/
theorem T_mulVec_D_xi (D T : Matrix n n ℝ) (β η ξ : n → ℝ)
    (hcomm : D * T - T * D = Matrix.vecMulVec β η - Matrix.vecMulVec η β)
    (hTxi : T *ᵥ ξ = 0) (hη : η ⬝ᵥ ξ = 1) (hβ : β ⬝ᵥ ξ = 0) :
    T *ᵥ (D *ᵥ ξ) = -β := by
  have rhs : (D * T - T * D) *ᵥ ξ = β := by
    rw [hcomm, Matrix.sub_mulVec, vecMulVec_act, vecMulVec_act, hη, hβ, one_smul, zero_smul,
      sub_zero]
  have lhs : (D * T - T * D) *ᵥ ξ = - (T *ᵥ (D *ᵥ ξ)) := by
    rw [Matrix.sub_mulVec, ← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, hTxi,
      Matrix.mulVec_zero, zero_sub]
  have h : -(T *ᵥ (D *ᵥ ξ)) = β := lhs.symm.trans rhs
  exact neg_eq_iff_eq_neg.mp h

#print axioms T_mulVec_D_xi

end CMKeyLemma
