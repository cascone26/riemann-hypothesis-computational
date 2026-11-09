import CCMConfluentPerp

/-!
# S2994 — Schur test for entrywise-perturbed quadratic forms (lane 2: keeps per-pair |ν_i − ν_k|)

* `schur_quad` — |M_ik − P_ik| ≤ δ_ik, δ symmetric, every row sum Σ_k δ_ik ≤ S ⇒ |xᵀMx − xᵀPx| ≤ S·Σx².
Applied with P = u𝟙ᵀ + 𝟙uᵀ (vanishes on 𝟙^⊥, S2990) and δ_ik = D₂|ν_i − ν_k|/π (S2990 entrywise), the 𝟙^⊥ tail bound becomes
(D₂/π)·max_i Σ_k|ν_i − ν_k|·‖x‖² instead of (2V D₂/π)·n·‖x‖² (S2992): a factor ≈ 2 for equispaced ν.
-/

open Finset Matrix

namespace CCMSchur

variable {n : Type*} [Fintype n]

theorem schur_quad (M P δ : Matrix n n ℝ) (S : ℝ) (h : ∀ i k, |M i k - P i k| ≤ δ i k)
    (hsym : ∀ i k, δ i k = δ k i) (hS : ∀ i, ∑ k, δ i k ≤ S) (x : n → ℝ) :
    |x ⬝ᵥ (M *ᵥ x) - x ⬝ᵥ (P *ᵥ x)| ≤ S * ∑ i, x i ^ 2 := by
  have hδ0 : ∀ i k, 0 ≤ δ i k := fun i k => le_trans (abs_nonneg _) (h i k)
  have e1 : x ⬝ᵥ (M *ᵥ x) - x ⬝ᵥ (P *ᵥ x) = ∑ i, ∑ k, x i * (M i k - P i k) * x k := by
    simp only [dotProduct, mulVec, Finset.mul_sum, ← Finset.sum_sub_distrib]
    apply Finset.sum_congr rfl; intro i _
    apply Finset.sum_congr rfl; intro k _
    ring
  rw [e1]
  have step1 : |∑ i, ∑ k, x i * (M i k - P i k) * x k| ≤ ∑ i, ∑ k, δ i k * ((x i ^ 2 + x k ^ 2) / 2) := by
    refine le_trans (Finset.abs_sum_le_sum_abs _ _) (Finset.sum_le_sum (fun i _ => ?_))
    refine le_trans (Finset.abs_sum_le_sum_abs _ _) (Finset.sum_le_sum (fun k _ => ?_))
    rw [abs_mul, abs_mul]
    have hxy : |x i| * |x k| ≤ (x i ^ 2 + x k ^ 2) / 2 := by
      nlinarith [sq_nonneg (|x i| - |x k|), sq_abs (x i), sq_abs (x k)]
    have := h i k
    calc |x i| * |M i k - P i k| * |x k| = |M i k - P i k| * (|x i| * |x k|) := by ring
      _ ≤ δ i k * ((x i ^ 2 + x k ^ 2) / 2) :=
          mul_le_mul this hxy (by positivity) (hδ0 i k)
  have split : ∑ i, ∑ k, δ i k * ((x i ^ 2 + x k ^ 2) / 2)
      = (∑ i, x i ^ 2 * ∑ k, δ i k) / 2 + (∑ k, x k ^ 2 * ∑ i, δ i k) / 2 := by
    rw [Finset.sum_comm (f := fun i k => δ i k * ((x i ^ 2 + x k ^ 2) / 2))]
    rw [show (∑ k, x k ^ 2 * ∑ i, δ i k) = ∑ k, ∑ i, x k ^ 2 * δ i k by
      apply Finset.sum_congr rfl; intro k _; rw [Finset.mul_sum]]
    rw [show (∑ i, x i ^ 2 * ∑ k, δ i k) = ∑ i, ∑ k, x i ^ 2 * δ i k by
      apply Finset.sum_congr rfl; intro i _; rw [Finset.mul_sum]]
    rw [Finset.sum_comm (f := fun i k => x i ^ 2 * δ i k)]
    rw [Finset.sum_div, Finset.sum_div, ← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl; intro k _
    rw [Finset.sum_div, Finset.sum_div, ← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl; intro i _
    ring
  have hrow : ∀ i, x i ^ 2 * ∑ k, δ i k ≤ x i ^ 2 * S :=
    fun i => mul_le_mul_of_nonneg_left (hS i) (sq_nonneg _)
  have hcol : ∀ k, x k ^ 2 * ∑ i, δ i k ≤ x k ^ 2 * S := by
    intro k
    have : ∑ i, δ i k = ∑ i, δ k i := Finset.sum_congr rfl (fun i _ => hsym i k)
    rw [this]; exact mul_le_mul_of_nonneg_left (hS k) (sq_nonneg _)
  have h1 := Finset.sum_le_sum (fun i (_ : i ∈ univ) => hrow i)
  have h2 := Finset.sum_le_sum (fun k (_ : k ∈ univ) => hcol k)
  rw [← Finset.sum_mul] at h1 h2
  calc |∑ i, ∑ k, x i * (M i k - P i k) * x k|
      ≤ ∑ i, ∑ k, δ i k * ((x i ^ 2 + x k ^ 2) / 2) := step1
    _ = (∑ i, x i ^ 2 * ∑ k, δ i k) / 2 + (∑ k, x k ^ 2 * ∑ i, δ i k) / 2 := split
    _ ≤ (∑ i, x i ^ 2) * S / 2 + (∑ k, x k ^ 2) * S / 2 := by linarith
    _ = S * ∑ i, x i ^ 2 := by ring

end CCMSchur

#print axioms CCMSchur.schur_quad
