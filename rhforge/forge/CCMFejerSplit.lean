import CCMWeightShape
import CCMLoewner
import CCMTruncMono

/-!
# S2951 — linearity of the zero-side Gram in the weights + Fejér split of the CCM zero side (lane 2)

* `gram_add`, `gram_smul` : W ↦ C·diag(W)·Cᵀ is linear;
* `gram_trunc_split` : G_T = G_S + G_{T∖S}  (S ⊆ T; truncated weights split additively);
* `ccm_gram_fejer` : G_ccm = (2/L)·C·Cᵀ − (2/L)·C·diag(cos γ_j L)·Cᵀ — from S2943 w_γ = 2(1−cos γL)/L:
  the zero side = a NON-oscillating Cauchy Gram (2/L)·CCᵀ minus an OSCILLATORY Gram carrying all the γ-dependence.
-/

open Matrix Finset

namespace CCMFejer

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

theorem gram_add (C : Matrix n m ℝ) (w w' : m → ℝ) :
    C * diagonal (w + w') * Cᵀ = C * diagonal w * Cᵀ + C * diagonal w' * Cᵀ := by
  rw [show diagonal (w + w') = diagonal w + diagonal w' from (diagonal_add w w').symm,
    Matrix.mul_add, Matrix.add_mul]

theorem gram_smul (C : Matrix n m ℝ) (a : ℝ) (w : m → ℝ) :
    C * diagonal (a • w) * Cᵀ = a • (C * diagonal w * Cᵀ) := by
  rw [show diagonal (a • w) = a • diagonal w by ext i j; by_cases h : i = j <;> simp [diagonal_apply, h]]
  rw [Matrix.mul_smul, Matrix.smul_mul]

theorem truncW_split {S T : Finset m} (hST : S ⊆ T) (w : m → ℝ) :
    CCMTrunc.truncW T w = CCMTrunc.truncW S w + CCMTrunc.truncW (T \ S) w := by
  funext j
  unfold CCMTrunc.truncW
  by_cases hS : j ∈ S
  · simp [hS, hST hS]
  · by_cases hT : j ∈ T <;> simp [hS, hT]

theorem gram_trunc_split (C : Matrix n m ℝ) {S T : Finset m} (hST : S ⊆ T) (w : m → ℝ) :
    C * diagonal (CCMTrunc.truncW T w) * Cᵀ
      = C * diagonal (CCMTrunc.truncW S w) * Cᵀ + C * diagonal (CCMTrunc.truncW (T \ S) w) * Cᵀ := by
  rw [truncW_split hST, gram_add]

theorem ccm_gram_fejer {L : ℝ} (C : Matrix n m ℝ) (γ : m → ℝ) :
    C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ
      = (2 / L) • (C * Cᵀ) - (2 / L) • (C * diagonal (fun j => Real.cos (γ j * L)) * Cᵀ) := by
  have hw : (fun j => CCMWeight.ccmW L (γ j))
      = (fun _ : m => 2 / L) + (-(2 / L)) • (fun j => Real.cos (γ j * L)) := by
    funext j
    simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, CCMWeight.ccmW_eq_cos]
    ring
  rw [hw, gram_add, gram_smul, CCMLoewner.const_gram, neg_smul]
  abel

end CCMFejer

#print axioms CCMFejer.gram_trunc_split
#print axioms CCMFejer.ccm_gram_fejer
