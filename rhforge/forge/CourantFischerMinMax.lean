/-
CourantFischerMinMax.lean (HP RH center, LEDGER S2898, 2026-11-07)
Courant–Fischer half 2 (dimension counting), building on RayleighEigenTail (S2897). Mathlib lacks this.
For symmetric T (eigenvalues sorted decreasingly):
  dim W ≥ k+1  ⟹ ∃ x ∈ W, x ≠ 0, ⟪Tx,x⟫ ≤ λ_k ‖x‖²
  dim W ≥ n−k  ⟹ ∃ x ∈ W, x ≠ 0, λ_k ‖x‖² ≤ ⟪Tx,x⟫
These are the two inequalities behind Cauchy interlacing (S2884's one-step interlacing).
-/
import Mathlib
import RayleighEigenTail

open Module

namespace CourantFischer

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [FiniteDimensional ℝ E]

local notation "⟪" x ", " y "⟫" => @inner ℝ E _ x y

/-- A subspace of dimension larger than |S| contains a nonzero vector orthogonal to all `b i`, `i ∈ S`. -/
theorem exists_orth_of_card_lt {n : ℕ} (b : Fin n → E) (S : Finset (Fin n)) (W : Submodule ℝ E)
    (h : S.card < finrank ℝ W) : ∃ x ∈ W, x ≠ 0 ∧ ∀ i ∈ S, ⟪b i, x⟫ = 0 := by
  let f : W →ₗ[ℝ] (S → ℝ) := LinearMap.pi fun j => (innerₛₗ ℝ (b j.1)).comp W.subtype
  have hrn := LinearMap.finrank_range_add_finrank_ker f
  have hrange : finrank ℝ (LinearMap.range f) ≤ S.card := by
    calc finrank ℝ (LinearMap.range f) ≤ finrank ℝ (S → ℝ) := Submodule.finrank_le _
      _ = S.card := by simp
  have hker : LinearMap.ker f ≠ ⊥ := by
    intro hbot
    rw [hbot, finrank_bot] at hrn
    omega
  obtain ⟨y, hy, hy0⟩ := Submodule.exists_mem_ne_zero_of_ne_bot hker
  refine ⟨(y : E), y.2, fun h0 => hy0 (Subtype.ext (by simpa using h0)), fun i hi => ?_⟩
  have := congrFun (LinearMap.mem_ker.mp hy) ⟨i, hi⟩
  simpa [f] using this

variable {n : ℕ} {T : E →ₗ[ℝ] E}

theorem exists_le_eigenvalue (hT : T.IsSymmetric) (hn : finrank ℝ E = n) (k : Fin n)
    (W : Submodule ℝ E) (hW : (k : ℕ) + 1 ≤ finrank ℝ W) :
    ∃ x ∈ W, x ≠ 0 ∧ ⟪T x, x⟫ ≤ hT.eigenvalues hn k * ‖x‖ ^ 2 := by
  obtain ⟨x, hxW, hx0, horth⟩ := exists_orth_of_card_lt (hT.eigenvectorBasis hn) (Finset.Iio k) W
    (by rw [Fin.card_Iio]; omega)
  exact ⟨x, hxW, hx0, RayleighTail.rayleigh_tail_le hT hn k x
    (fun i hi => horth i (Finset.mem_Iio.mpr hi))⟩

theorem exists_ge_eigenvalue (hT : T.IsSymmetric) (hn : finrank ℝ E = n) (k : Fin n)
    (W : Submodule ℝ E) (hW : n - k ≤ finrank ℝ W) :
    ∃ x ∈ W, x ≠ 0 ∧ hT.eigenvalues hn k * ‖x‖ ^ 2 ≤ ⟪T x, x⟫ := by
  obtain ⟨x, hxW, hx0, horth⟩ := exists_orth_of_card_lt (hT.eigenvectorBasis hn) (Finset.Ioi k) W
    (by rw [Fin.card_Ioi]; have := k.2; omega)
  exact ⟨x, hxW, hx0, RayleighTail.rayleigh_head_ge hT hn k x
    (fun i hi => horth i (Finset.mem_Ioi.mpr hi))⟩

end CourantFischer

#print axioms CourantFischer.exists_orth_of_card_lt
#print axioms CourantFischer.exists_le_eigenvalue
#print axioms CourantFischer.exists_ge_eigenvalue
