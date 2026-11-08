import CCMWeightGram

/-!
# S2937 — Loewner-order brackets for the CCM zero-side Gram (lane 2)

The weighted Cauchy Gram G_w = C·diag(w)·Cᵀ is MONOTONE in the weight vector in the Loewner (PSD) order:
w ≤ w' pointwise ⇒ G_{w'} − G_w = C·diag(w'−w)·Cᵀ ⪰ 0.  Specialised to the CCM/Shi weights
w_γ = 4 sin²(γL/2)/L ∈ [0, 4/L], this gives the matrix sandwich
  0 ⪯ G_ccm ⪯ (4/L)·C·Cᵀ,
upgrading S2931's scalar form bound (`ccm_form_le`) to a full operator inequality: every eigenvalue and Rayleigh
quantity of the zero side is dominated by that of the unweighted Cauchy Gram scaled by 4/L.
-/

open Matrix Finset

namespace CCMLoewner

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

theorem gram_sub (C : Matrix n m ℝ) (w w' : m → ℝ) :
    C * diagonal w' * Cᵀ - C * diagonal w * Cᵀ = C * diagonal (w' - w) * Cᵀ := by
  rw [← Matrix.sub_mul, ← Matrix.mul_sub, diagonal_sub]; rfl

theorem gram_mono (C : Matrix n m ℝ) {w w' : m → ℝ} (h : ∀ j, w j ≤ w' j) :
    (C * diagonal w' * Cᵀ - C * diagonal w * Cᵀ).PosSemidef := by
  rw [gram_sub]
  exact GramWeil.posSemidef_of_nonneg C _ (fun j => by simp [sub_nonneg.mpr (h j)])

theorem ccm_gram_loewner_le {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ) :
    (C * diagonal (fun _ : m => 4 / L) * Cᵀ - C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ).PosSemidef :=
  gram_mono C (fun j => CCMWeight.ccmW_le hL (γ j))

theorem const_gram (C : Matrix n m ℝ) (a : ℝ) :
    C * diagonal (fun _ : m => a) * Cᵀ = a • (C * Cᵀ) := by
  rw [← Matrix.smul_one_eq_diagonal]; simp [Matrix.mul_smul, Matrix.smul_mul]

/-- The sandwich 0 ⪯ G_ccm ⪯ (4/L)·C·Cᵀ. -/
theorem ccm_gram_sandwich {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ) :
    (C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ).PosSemidef ∧
    ((4 / L) • (C * Cᵀ) - C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ).PosSemidef := by
  refine ⟨CCMWeight.ccm_gram_posSemidef hL C γ, ?_⟩
  rw [← const_gram]
  exact ccm_gram_loewner_le hL C γ

end CCMLoewner

#print axioms CCMLoewner.gram_mono
#print axioms CCMLoewner.ccm_gram_sandwich
