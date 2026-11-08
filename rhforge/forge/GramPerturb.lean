/-
GramPerturb.lean (HP RH center, LEDGER S2927, 2026-11-08) — lane 2.
Robustness of the zero-side Weil form under weight perturbation: for G(w) = C·diag(w)·Cᵀ,
  |xᵀG(w')x − xᵀG(w)x| ≤ (Σ_j δ_j Σ_i C_ij²)·‖x‖²   whenever |w'_j − w_j| ≤ δ_j.
Consequence (`margin_stable`): a certified lower bound μ on the form at w transfers to w' with loss
Δ = Σ_j δ_j Σ_i C_ij² — the Lean form of "numerical error in the zero data moves the margin by a bounded amount".
-/
import Mathlib
import GramTruncation

open Matrix Finset

namespace GramPert

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq m]

open GramWeil

theorem form_sub (C : Matrix n m ℝ) (w w' : m → ℝ) (x : n → ℝ) :
    x ⬝ᵥ ((C * diagonal w' * Cᵀ) *ᵥ x) - x ⬝ᵥ ((C * diagonal w * Cᵀ) *ᵥ x)
      = ∑ j, (w' j - w j) * ((Cᵀ *ᵥ x) j) ^ 2 := by
  rw [quad_form, quad_form, ← sum_sub_distrib]
  exact sum_congr rfl fun j _ => by ring

theorem form_perturb (C : Matrix n m ℝ) (w w' δ : m → ℝ) (h : ∀ j, |w' j - w j| ≤ δ j)
    (x : n → ℝ) :
    |x ⬝ᵥ ((C * diagonal w' * Cᵀ) *ᵥ x) - x ⬝ᵥ ((C * diagonal w * Cᵀ) *ᵥ x)|
      ≤ (∑ j, δ j * ∑ i, C i j ^ 2) * ∑ i, x i ^ 2 := by
  rw [form_sub, Finset.sum_mul]
  refine (abs_sum_le_sum_abs _ _).trans (sum_le_sum fun j _ => ?_)
  rw [abs_mul, abs_of_nonneg (sq_nonneg ((Cᵀ *ᵥ x) j)), mul_assoc]
  have hδ : 0 ≤ δ j := (abs_nonneg _).trans (h j)
  exact mul_le_mul (h j) (GramTrunc.col_sq_le C x j) (sq_nonneg _) hδ

/-- **Margin stability.** -/
theorem margin_stable (C : Matrix n m ℝ) (w w' δ : m → ℝ) (h : ∀ j, |w' j - w j| ≤ δ j)
    (μ : ℝ) (hμ : ∀ x : n → ℝ, μ * ∑ i, x i ^ 2 ≤ x ⬝ᵥ ((C * diagonal w * Cᵀ) *ᵥ x))
    (x : n → ℝ) :
    (μ - ∑ j, δ j * ∑ i, C i j ^ 2) * ∑ i, x i ^ 2 ≤ x ⬝ᵥ ((C * diagonal w' * Cᵀ) *ᵥ x) := by
  have h1 := hμ x
  have h2 := (abs_le.mp (form_perturb C w w' δ h x)).1
  nlinarith [h1, h2]

end GramPert

#print axioms GramPert.form_sub
#print axioms GramPert.form_perturb
#print axioms GramPert.margin_stable
