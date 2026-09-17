import Mathlib
import EffectiveRank

/-!
# Matrix-level effective-rank bound (RH session 2026-08-29)

Lift of `effective_rank_bound` to a Hermitian real matrix: `(tr A)² ≤ rank A · ∑ λᵢ²`,
i.e. `rank A ≥ (tr A)² / ‖A‖²_HS`.  Uses `IsHermitian.trace_eq_sum_eigenvalues` +
`rank_eq_card_non_zero_eigs`.  Companion to `rank_trace_matrix`; completes the matrix toolkit.
-/

open Matrix Finset

/-- `(tr A)² ≤ rank A · ∑ λᵢ²` for a Hermitian real matrix. -/
theorem effective_rank_matrix {n : Type*} [Fintype n] [DecidableEq n] {A : Matrix n n ℝ}
    (hA : A.IsHermitian) :
    (A.trace) ^ 2 ≤ (A.rank : ℝ) * (∑ i, (hA.eigenvalues i) ^ 2) := by
  have ht : A.trace = ∑ i, hA.eigenvalues i := by
    simpa using hA.trace_eq_sum_eigenvalues
  have hr : (A.rank : ℝ)
      = ((univ.filter (fun i => hA.eigenvalues i ≠ 0)).card : ℝ) := by
    rw [hA.rank_eq_card_non_zero_eigs, Fintype.card_subtype]
  rw [ht, hr]
  exact effective_rank_bound univ hA.eigenvalues
