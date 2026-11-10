import CCMMarginSchur
import GramWeilForm

/-!
# S3005 — the far tail must be used EXACTLY, not bounded (lane 2: negative scan → reduction)

NEGATIVE (HANDLE numerics, S3005 scan): for λ² ∈ {3,5,10,20}, N ≤ 6 the near block never has ≥ n−1 zeros below 2V, and the
best tail bound (S3001, ~10⁻³–10⁻¹) exceeds the true margin e^{−12.3λ²} (10⁻¹⁷–10⁻¹⁰⁷) by 15–100 orders.  Bounding the tail
cannot certify CCM margins.  The correct use: the tail is a Gram matrix with weights w_γ = ccmW ≥ 0 exactly when the far
zeros are real, so it is PSD and can only HELP:
* `tail_gram_nonneg` — w ≥ 0 ⇒ xᵀ(C diag(w) Cᵀ)x ≥ 0.
* `margin_from_near_and_real_tail` — xᵀAx = xᵀG_Kx + xᵀ(C_f diag(w_f) C_fᵀ)x, w_f ≥ 0 ⇒ xᵀG_Kx ≤ xᵀAx.
Reading (Report): CCM positivity of A on a window ⇐ (finite verified zeros give G_K ⪰ m) + (zeros above T real ⇒ tail PSD);
the second clause is exactly RH above height T — the reduction is sharp and makes clear where the content of RH sits.
-/

open Matrix Finset

namespace CCMTailExact

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq m]

theorem tail_gram_nonneg (C : Matrix n m ℝ) (w : m → ℝ) (hw : ∀ j, 0 ≤ w j) (x : n → ℝ) :
    0 ≤ x ⬝ᵥ ((C * diagonal w * Cᵀ) *ᵥ x) := by
  have h := (GramWeil.posSemidef_of_nonneg C w hw).dotProduct_mulVec_nonneg x
  simpa using h

theorem margin_from_near_and_real_tail (A G : Matrix n n ℝ) (C : Matrix n m ℝ) (w : m → ℝ)
    (hw : ∀ j, 0 ≤ w j) (x : n → ℝ)
    (hsplit : x ⬝ᵥ (A *ᵥ x) = x ⬝ᵥ (G *ᵥ x) + x ⬝ᵥ ((C * diagonal w * Cᵀ) *ᵥ x)) :
    x ⬝ᵥ (G *ᵥ x) ≤ x ⬝ᵥ (A *ᵥ x) := by
  have := tail_gram_nonneg C w hw x
  linarith

end CCMTailExact

#print axioms CCMTailExact.tail_gram_nonneg
#print axioms CCMTailExact.margin_from_near_and_real_tail
