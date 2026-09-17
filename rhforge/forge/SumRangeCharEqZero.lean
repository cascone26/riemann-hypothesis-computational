import Mathlib.NumberTheory.LSeries.Dirichlet
import Mathlib.Data.ZMod.Basic

open Finset

theorem sum_range_char_eq_zero {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M} (hχ : χ ≠ 1) :
    ∑ k ∈ Finset.range M, χ (k : ZMod M) = 0 := by
  rw [← MulChar.sum_eq_zero_of_ne_one hχ]
  refine Finset.sum_nbij' (fun k => (k : ZMod M)) (fun a => a.val) ?_ ?_ ?_ ?_ ?_
  · intro k _; exact Finset.mem_univ _
  · intro a _; exact Finset.mem_range.mpr (ZMod.val_lt a)
  · intro k hk; exact ZMod.val_cast_of_lt (Finset.mem_range.mp hk)
  · intro a _; exact ZMod.natCast_rightInverse a
  · intro k _; rfl

#print axioms sum_range_char_eq_zero
