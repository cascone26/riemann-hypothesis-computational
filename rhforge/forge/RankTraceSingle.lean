import Mathlib

/-!
# Single-matrix rank–trace shadow (RH session 2026-08-29)

The `P`-part of the Alpöge–Furman rank–trace inequality (AF eq (L); scalar integrality
`(p-1)^2 ≥ 0`, line 587):  for any real spectrum `p`,
`rank P = #{i : pᵢ ≠ 0} ≥ 2·tr P − ‖P‖²_HS = 2 ∑ pᵢ − ∑ pᵢ²`.
Per eigenvalue `[pᵢ≠0] ≥ 2pᵢ − pᵢ²`, since `1 - (2p - p²) = (p-1)² ≥ 0` (and `0 ≥ 0` at `p=0`).
-/

open Finset

/-- `2 ∑ pᵢ − ∑ pᵢ² ≤ #{i : pᵢ ≠ 0}` for any real family `p`. -/
theorem rank_trace_single {ι : Type*} (s : Finset ι) (p : ι → ℝ) :
    2 * (∑ i ∈ s, p i) - (∑ i ∈ s, p i ^ 2)
      ≤ (s.filter (fun i => p i ≠ 0)).card := by
  classical
  have per : ∀ i ∈ s, 2 * p i - p i ^ 2 ≤ (if p i ≠ 0 then (1:ℝ) else 0) := by
    intro i _
    by_cases h : p i = 0
    · simp [h]
    · rw [if_pos h]; nlinarith [sq_nonneg (p i - 1)]
  have hsum : ∑ i ∈ s, (2 * p i - p i ^ 2)
      ≤ ∑ i ∈ s, (if p i ≠ 0 then (1:ℝ) else 0) := Finset.sum_le_sum per
  have hL : ∑ i ∈ s, (2 * p i - p i ^ 2)
      = 2 * (∑ i ∈ s, p i) - (∑ i ∈ s, p i ^ 2) := by
    simp only [Finset.sum_sub_distrib, Finset.mul_sum]
  have hR : ∑ i ∈ s, (if p i ≠ 0 then (1:ℝ) else 0)
      = ((s.filter (fun i => p i ≠ 0)).card : ℝ) := by
    rw [Finset.sum_ite, Finset.sum_const, Finset.sum_const_zero, add_zero,
        nsmul_eq_mul, mul_one]
  rw [hL, hR] at hsum
  exact hsum
