import Mathlib.NumberTheory.LSeries.Dirichlet
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv

open Finset

/-! §1525 the twisted partial sums are entire (hF input for the holomorphic-limit step). EXIT0 / std-axioms. -/

theorem differentiable_twist_partial_sum {N : ℕ} [NeZero N] (χ : DirichletCharacter ℂ N) (n : ℕ) :
    Differentiable ℂ (fun s : ℂ => ∑ k ∈ Finset.Icc 1 n, χ (k : ZMod N) * (k : ℂ) ^ (-s)) := by
  have heq : (fun s : ℂ => ∑ k ∈ Finset.Icc 1 n, χ (k : ZMod N) * (k : ℂ) ^ (-s))
      = ∑ k ∈ Finset.Icc 1 n, (fun s : ℂ => χ (k : ZMod N) * (k : ℂ) ^ (-s)) := by
    funext s; rw [Finset.sum_apply]
  rw [heq]
  apply Differentiable.sum
  intro k hk
  rw [Finset.mem_Icc] at hk
  have hk0 : (k : ℂ) ≠ 0 := by
    have h1 : (1 : ℕ) ≤ k := hk.1
    exact_mod_cast Nat.one_le_iff_ne_zero.mp h1
  exact (differentiable_const _).mul (Differentiable.neg differentiable_id |>.const_cpow (Or.inl hk0))

#print axioms differentiable_twist_partial_sum
