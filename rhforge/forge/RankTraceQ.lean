import Mathlib

/-!
# Q-part rank–trace shadow (RH session 2026-08-29)

The `Q`-part of Alpöge–Furman's rank–trace inequality (AF eq (L); the "flat charge 4" integrality
`(q-2)^2 ≥ 0`, line 587):  for any real spectrum `q`,
`4 ∑ qᵢ − ∑ qᵢ² ≤ 4 · #{i : qᵢ > 0} ( ≤ 4·npl(Q) )`.
Per eigenvalue `4qᵢ − qᵢ² ≤ 4·[qᵢ>0]`: for `q>0`, `4 - (4q - q²) = (q-2)² ≥ 0`; for `q ≤ 0`,
`4q - q² ≤ 0`.  Pairs with `RankTraceSingle` toward the full `r ≥ 2trP + 4trQ − 4b − ‖P+Q‖²`.
-/

open Finset

/-- `4 ∑ qᵢ − ∑ qᵢ² ≤ 4 · #{i : qᵢ > 0}` for any real family `q`. -/
theorem rank_trace_Q {ι : Type*} (s : Finset ι) (q : ι → ℝ) :
    4 * (∑ i ∈ s, q i) - (∑ i ∈ s, q i ^ 2)
      ≤ 4 * (s.filter (fun i => 0 < q i)).card := by
  classical
  have per : ∀ i ∈ s, 4 * q i - q i ^ 2 ≤ (if 0 < q i then (4:ℝ) else 0) := by
    intro i _
    by_cases h : 0 < q i
    · rw [if_pos h]; nlinarith [sq_nonneg (q i - 2)]
    · rw [if_neg h]; have hle : q i ≤ 0 := not_lt.mp h
      nlinarith [sq_nonneg (q i), hle]
  have hsum : ∑ i ∈ s, (4 * q i - q i ^ 2)
      ≤ ∑ i ∈ s, (if 0 < q i then (4:ℝ) else 0) := Finset.sum_le_sum per
  have hL : ∑ i ∈ s, (4 * q i - q i ^ 2)
      = 4 * (∑ i ∈ s, q i) - (∑ i ∈ s, q i ^ 2) := by
    simp only [Finset.sum_sub_distrib, Finset.mul_sum]
  have hR : ∑ i ∈ s, (if 0 < q i then (4:ℝ) else 0)
      = 4 * ((s.filter (fun i => 0 < q i)).card : ℝ) := by
    rw [Finset.sum_ite, Finset.sum_const, Finset.sum_const_zero, add_zero,
        nsmul_eq_mul, mul_comm]
  rw [hL, hR] at hsum
  exact hsum
