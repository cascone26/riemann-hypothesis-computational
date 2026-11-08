/-
RayleighEigenTail.lean (HP RH center, LEDGER S2897, 2026-11-07)
First half of Courant–Fischer (absent from Mathlib), toward a full Lean proof of S2884's one-step
interlacing. For a symmetric T on a finite-dimensional real inner product space, with Mathlib's
eigenvalues sorted DECREASINGLY (eigenvalues_antitone) and orthonormal eigenvectorBasis b:
  (tail)  x ⊥ b i for all i < k   ⟹  ⟪Tx,x⟫ ≤ λ_k ‖x‖²
  (head)  x ⊥ b i for all i > k   ⟹  λ_k ‖x‖² ≤ ⟪Tx,x⟫
-/
import Mathlib

open Module

namespace RayleighTail

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [FiniteDimensional ℝ E]
  {n : ℕ} {T : E →ₗ[ℝ] E}

local notation "⟪" x ", " y "⟫" => @inner ℝ E _ x y

lemma inner_T_eq_sum (hT : T.IsSymmetric) (hn : finrank ℝ E = n) (x : E) :
    ⟪T x, x⟫ = ∑ i, hT.eigenvalues hn i * ⟪hT.eigenvectorBasis hn i, x⟫ ^ 2 := by
  set b := hT.eigenvectorBasis hn
  nth_rewrite 1 [← b.sum_repr' x]
  rw [map_sum, sum_inner]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [map_smul, hT.apply_eigenvectorBasis, smul_smul, inner_smul_left]
  simp [b, sq, mul_comm, mul_left_comm]

lemma norm_sq_eq_sum (hT : T.IsSymmetric) (hn : finrank ℝ E = n) (x : E) :
    ‖x‖ ^ 2 = ∑ i, ⟪hT.eigenvectorBasis hn i, x⟫ ^ 2 := by
  rw [← (hT.eigenvectorBasis hn).sum_sq_norm_inner_right x]
  simp [Real.norm_eq_abs, sq_abs]

theorem rayleigh_tail_le (hT : T.IsSymmetric) (hn : finrank ℝ E = n) (k : Fin n) (x : E)
    (hx : ∀ i : Fin n, i < k → ⟪hT.eigenvectorBasis hn i, x⟫ = 0) :
    ⟪T x, x⟫ ≤ hT.eigenvalues hn k * ‖x‖ ^ 2 := by
  rw [inner_T_eq_sum hT hn, norm_sq_eq_sum hT hn, Finset.mul_sum]
  refine Finset.sum_le_sum fun i _ => ?_
  rcases lt_or_ge i k with h | h
  · simp [hx i h]
  · exact mul_le_mul_of_nonneg_right (hT.eigenvalues_antitone hn h) (sq_nonneg _)

theorem rayleigh_head_ge (hT : T.IsSymmetric) (hn : finrank ℝ E = n) (k : Fin n) (x : E)
    (hx : ∀ i : Fin n, k < i → ⟪hT.eigenvectorBasis hn i, x⟫ = 0) :
    hT.eigenvalues hn k * ‖x‖ ^ 2 ≤ ⟪T x, x⟫ := by
  rw [inner_T_eq_sum hT hn, norm_sq_eq_sum hT hn, Finset.mul_sum]
  refine Finset.sum_le_sum fun i _ => ?_
  rcases lt_or_ge k i with h | h
  · simp [hx i h]
  · exact mul_le_mul_of_nonneg_right (hT.eigenvalues_antitone hn h) (sq_nonneg _)

end RayleighTail

#print axioms RayleighTail.rayleigh_tail_le
#print axioms RayleighTail.rayleigh_head_ge
