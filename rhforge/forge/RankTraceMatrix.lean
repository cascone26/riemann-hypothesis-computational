import Mathlib
import RankTraceSingle

/-!
# Matrix-level rank–trace shadow (RH session 2026-08-29)

Lift of `rank_trace_single` to an actual Hermitian real matrix: for `A` Hermitian,
`2 · tr A − ∑ λᵢ² ≤ rank A`, where `λᵢ = A.eigenvalues` and `∑λᵢ² = ‖A‖²_HS`.
Uses Mathlib's `IsHermitian.trace_eq_sum_eigenvalues` and `rank_eq_card_non_zero_eigs`.
-/

open Matrix Finset

/-- Matrix form: `2 · tr A − ∑ λᵢ² ≤ rank A` for a Hermitian real matrix. -/
theorem rank_trace_matrix {n : Type*} [Fintype n] [DecidableEq n] {A : Matrix n n ℝ}
    (hA : A.IsHermitian) :
    2 * A.trace - (∑ i, (hA.eigenvalues i) ^ 2) ≤ (A.rank : ℝ) := by
  have ht : A.trace = ∑ i, hA.eigenvalues i := by
    simpa using hA.trace_eq_sum_eigenvalues
  have hr : (A.rank : ℝ)
      = ((univ.filter (fun i => hA.eigenvalues i ≠ 0)).card : ℝ) := by
    rw [hA.rank_eq_card_non_zero_eigs, Fintype.card_subtype]
  rw [ht, hr]
  exact rank_trace_single univ hA.eigenvalues
