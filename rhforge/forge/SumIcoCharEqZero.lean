import Mathlib.NumberTheory.LSeries.Dirichlet
import Mathlib.Data.ZMod.Basic

open Finset

/-!
# §1513 — shifted-block character sum: `∑_{k ∈ Ico n (n+M)} χ(↑k) = 0` (`χ ≠ 1`)

Any `M` consecutive integers hit each residue mod `M` exactly once, so a full block sums to the full period sum
`∑_{a : ZMod M} χ a = 0` (`MulChar.sum_eq_zero_of_ne_one`).  Reindex `Ico n (n+M) → range M` by the shift, then the
translated bijection `j ↦ ↑(n+j)`.  This is the periodicity input `S(n+M) = S(n)` for the partial sums feeding
Abel summation.  EXIT0 / std-axioms.
-/

theorem sum_Ico_char_eq_zero {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M} (hχ : χ ≠ 1) (n : ℕ) :
    ∑ k ∈ Finset.Ico n (n + M), χ (k : ZMod M) = 0 := by
  rw [Finset.sum_Ico_eq_sum_range, Nat.add_sub_cancel_left, ← MulChar.sum_eq_zero_of_ne_one hχ]
  refine Finset.sum_nbij' (fun j => ((n + j : ℕ) : ZMod M)) (fun a => (a - (n : ZMod M)).val)
    ?_ ?_ ?_ ?_ ?_
  · intro j _; exact Finset.mem_univ _
  · intro a _; exact Finset.mem_range.mpr (ZMod.val_lt _)
  · intro j hj
    have h1 : ((n + j : ℕ) : ZMod M) - (n : ZMod M) = (j : ZMod M) := by push_cast; ring
    rw [h1]; exact ZMod.val_cast_of_lt (Finset.mem_range.mp hj)
  · intro a _
    rw [Nat.cast_add, ZMod.natCast_rightInverse (a - (n : ZMod M))]; ring
  · intro j _; rfl

#print axioms sum_Ico_char_eq_zero
