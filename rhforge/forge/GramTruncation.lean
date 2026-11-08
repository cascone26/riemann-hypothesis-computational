/-
GramTruncation.lean (HP RH center, LEDGER S2925, 2026-11-08) — lane 2.
Controlling the truncation of the zero sum in the CCM/Shi zero-side form xᵀ(C·diag(w)·Cᵀ)x:
- `form_mono`: w ≤ w' entrywise ⟹ form(w) ≤ form(w')   (Loewner monotone in the weights);
- `form_add`: C·diag(w+w')·Cᵀ = C·diag(w)·Cᵀ + C·diag(w')·Cᵀ;
- `form_split`: form(w) = Σ_{j∈s} w_j(Cᵀx)_j² + Σ_{j∉s} w_j(Cᵀx)_j²   (kept zeros + tail);
- `tail_le`: w ≥ 0 ⟹ tail ≤ (Σ_{j∉s} w_j Σ_i C_ij²)·‖x‖²   (Cauchy–Schwarz).
So a computed partial zero sum plus an explicit column-norm tail bound brackets the full zero side:
form_s ≤ form ≤ form_s + T·‖x‖², with T = Σ_{tail} w_j Σ_i 1/(ν_i+γ_j)² = O(Σ_{γ>Γ} 1/γ²) in the dictionary.
-/
import Mathlib
import GramWeilForm

open Matrix Finset

namespace GramTrunc

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq m]

open GramWeil

theorem form_mono (C : Matrix n m ℝ) (w w' : m → ℝ) (h : ∀ j, w j ≤ w' j) (x : n → ℝ) :
    x ⬝ᵥ ((C * diagonal w * Cᵀ) *ᵥ x) ≤ x ⬝ᵥ ((C * diagonal w' * Cᵀ) *ᵥ x) := by
  rw [quad_form, quad_form]
  exact sum_le_sum fun j _ => mul_le_mul_of_nonneg_right (h j) (sq_nonneg _)

theorem form_add (C : Matrix n m ℝ) (w w' : m → ℝ) :
    C * diagonal (w + w') * Cᵀ = C * diagonal w * Cᵀ + C * diagonal w' * Cᵀ := by
  have : diagonal (w + w') = diagonal w + diagonal w' := (diagonal_add w w').symm
  rw [this, Matrix.mul_add, Matrix.add_mul]

theorem form_split (C : Matrix n m ℝ) (w : m → ℝ) (x : n → ℝ) (s : Finset m) :
    x ⬝ᵥ ((C * diagonal w * Cᵀ) *ᵥ x)
      = ∑ j ∈ s, w j * ((Cᵀ *ᵥ x) j) ^ 2 + ∑ j ∈ sᶜ, w j * ((Cᵀ *ᵥ x) j) ^ 2 := by
  rw [quad_form, sum_add_sum_compl]

theorem col_sq_le (C : Matrix n m ℝ) (x : n → ℝ) (j : m) :
    ((Cᵀ *ᵥ x) j) ^ 2 ≤ (∑ i, C i j ^ 2) * ∑ i, x i ^ 2 := by
  have : (Cᵀ *ᵥ x) j = ∑ i, C i j * x i := by
    simp [mulVec, dotProduct, transpose_apply]
  rw [this]
  exact Finset.sum_mul_sq_le_sq_mul_sq _ _ _

theorem tail_le (C : Matrix n m ℝ) (w : m → ℝ) (hw : ∀ j, 0 ≤ w j) (x : n → ℝ) (s : Finset m) :
    ∑ j ∈ sᶜ, w j * ((Cᵀ *ᵥ x) j) ^ 2
      ≤ (∑ j ∈ sᶜ, w j * ∑ i, C i j ^ 2) * ∑ i, x i ^ 2 := by
  rw [Finset.sum_mul]
  refine sum_le_sum fun j _ => ?_
  rw [mul_assoc]
  exact mul_le_mul_of_nonneg_left (col_sq_le C x j) (hw j)

/-- **Bracket**: kept part ≤ full form ≤ kept part + T‖x‖². -/
theorem bracket (C : Matrix n m ℝ) (w : m → ℝ) (hw : ∀ j, 0 ≤ w j) (x : n → ℝ) (s : Finset m) :
    ∑ j ∈ s, w j * ((Cᵀ *ᵥ x) j) ^ 2 ≤ x ⬝ᵥ ((C * diagonal w * Cᵀ) *ᵥ x) ∧
    x ⬝ᵥ ((C * diagonal w * Cᵀ) *ᵥ x)
      ≤ ∑ j ∈ s, w j * ((Cᵀ *ᵥ x) j) ^ 2 + (∑ j ∈ sᶜ, w j * ∑ i, C i j ^ 2) * ∑ i, x i ^ 2 := by
  rw [form_split C w x s]
  refine ⟨?_, ?_⟩
  · have : 0 ≤ ∑ j ∈ sᶜ, w j * ((Cᵀ *ᵥ x) j) ^ 2 :=
      sum_nonneg fun j _ => mul_nonneg (hw j) (sq_nonneg _)
    linarith
  · linarith [tail_le C w hw x s]

end GramTrunc

#print axioms GramTrunc.form_mono
#print axioms GramTrunc.form_add
#print axioms GramTrunc.tail_le
#print axioms GramTrunc.bracket
