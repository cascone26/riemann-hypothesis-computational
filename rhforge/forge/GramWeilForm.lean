/-
GramWeilForm.lean (HP RH center, LEDGER S2923, 2026-11-08) — lane 2 (Weil-positivity frontier).
Finite Weil-form of the CCM τ ≡ Shi S_{N,L} zero side. For G = C·diag(w)·Cᵀ (C any real n×m):
  xᵀ G x = Σ_j w_j (Cᵀx)_j²          (`quad_form`)
  w ≥ 0 ⟹ G positive semidefinite   (`psd_of_nonneg`)
  w > 0, Cᵀ injective ⟹ xᵀGx > 0 for x ≠ 0   (`pos_of_pos`)
  CERTIFICATE: ∃x, xᵀGx < 0 ⟹ ∃j, w_j < 0   (`neg_weight_of_neg_form`)
In the dictionary (C_mj = 1/(ν_m+γ_j), w_j = 4sin²(γ_j L/2)/L for on-line zeros), a single negative value of the
quadratic form of τ = W02 − W_R − W_p is a certificate that the weight structure fails — i.e. the finite
Weil/CCM criterion's contrapositive, stated as checkable linear algebra.
-/
import Mathlib

open Matrix Finset

namespace GramWeil

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq m]

theorem quad_form (C : Matrix n m ℝ) (w : m → ℝ) (x : n → ℝ) :
    x ⬝ᵥ ((C * diagonal w * Cᵀ) *ᵥ x) = ∑ j, w j * ((Cᵀ *ᵥ x) j) ^ 2 := by
  rw [← mulVec_mulVec, ← mulVec_mulVec, dotProduct_mulVec, ← mulVec_transpose]
  simp only [dotProduct, mulVec_diagonal]
  refine sum_congr rfl fun j _ => ?_
  ring

theorem psd_of_nonneg (C : Matrix n m ℝ) (w : m → ℝ) (hw : ∀ j, 0 ≤ w j) (x : n → ℝ) :
    0 ≤ x ⬝ᵥ ((C * diagonal w * Cᵀ) *ᵥ x) := by
  rw [quad_form]
  exact sum_nonneg fun j _ => mul_nonneg (hw j) (sq_nonneg _)

theorem pos_of_pos (C : Matrix n m ℝ) (w : m → ℝ) (hw : ∀ j, 0 < w j)
    (hinj : Function.Injective (Cᵀ).mulVec) (x : n → ℝ) (hx : x ≠ 0) :
    0 < x ⬝ᵥ ((C * diagonal w * Cᵀ) *ᵥ x) := by
  rw [quad_form]
  have hne : Cᵀ *ᵥ x ≠ 0 := by
    intro h; apply hx; apply hinj; rw [h, mulVec_zero]
  obtain ⟨j, hj⟩ : ∃ j, (Cᵀ *ᵥ x) j ≠ 0 := by
    by_contra h; push Not at h; exact hne (funext h)
  have hpos : 0 < w j * ((Cᵀ *ᵥ x) j) ^ 2 := mul_pos (hw j) (by positivity)
  exact lt_of_lt_of_le hpos (single_le_sum (f := fun j => w j * ((Cᵀ *ᵥ x) j) ^ 2)
    (fun i _ => mul_nonneg (hw i).le (sq_nonneg _)) (mem_univ j))

/-- **Certificate form**: a negative value of the quadratic form forces a negative weight. -/
theorem neg_weight_of_neg_form (C : Matrix n m ℝ) (w : m → ℝ) (x : n → ℝ)
    (hneg : x ⬝ᵥ ((C * diagonal w * Cᵀ) *ᵥ x) < 0) : ∃ j, w j < 0 := by
  by_contra h; push Not at h
  exact absurd (psd_of_nonneg C w h x) (not_le.mpr hneg)

/-- Mathlib `PosSemidef` packaging. -/
theorem posSemidef_of_nonneg (C : Matrix n m ℝ) (w : m → ℝ) (hw : ∀ j, 0 ≤ w j) :
    (C * diagonal w * Cᵀ).PosSemidef := by
  have hd : (diagonal w).PosSemidef := by
    rw [posSemidef_diagonal_iff]; intro j; exact_mod_cast hw j
  have := hd.mul_mul_conjTranspose_same C
  simpa [conjTranspose_eq_transpose_of_trivial] using this

end GramWeil

#print axioms GramWeil.quad_form
#print axioms GramWeil.psd_of_nonneg
#print axioms GramWeil.pos_of_pos
#print axioms GramWeil.neg_weight_of_neg_form
#print axioms GramWeil.posSemidef_of_nonneg
