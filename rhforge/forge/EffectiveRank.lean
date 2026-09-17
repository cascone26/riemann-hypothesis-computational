import Mathlib

/-!
# Effective-rank / Cauchy–Schwarz bound (RH session 2026-08-29, EP-Weil-discriminator work)

The engine behind the Alpöge–Furman rank–trace inequality: for a family `p : ι → ℝ`
(the eigenvalues of a PSD matrix P), the number of nonzero entries (rank P) satisfies
`(∑ pᵢ)² ≤ (#support) · (∑ pᵢ²)`, i.e. `rank P ≥ (tr P)² / ‖P‖²_HS`.  Pure Cauchy–Schwarz.
-/

open Finset

/-- Cauchy–Schwarz effective-rank bound: `(∑ pᵢ)² ≤ (#support) · ∑ pᵢ²`. -/
theorem effective_rank_bound {ι : Type*} (s : Finset ι) (p : ι → ℝ) :
    (∑ i ∈ s, p i) ^ 2
      ≤ (s.filter (fun i => p i ≠ 0)).card * (∑ i ∈ s, p i ^ 2) := by
  classical
  set t := s.filter (fun i => p i ≠ 0) with ht
  have hts : t ⊆ s := Finset.filter_subset _ _
  have hsum : ∑ i ∈ t, p i = ∑ i ∈ s, p i := by
    apply Finset.sum_subset hts
    intro i hi hnit
    by_contra hpne
    exact hnit (Finset.mem_filter.mpr ⟨hi, hpne⟩)
  have hsq : ∑ i ∈ t, p i ^ 2 ≤ ∑ i ∈ s, p i ^ 2 :=
    Finset.sum_le_sum_of_subset_of_nonneg hts (fun i _ _ => by positivity)
  have hcs : (∑ i ∈ t, (1:ℝ) * p i) ^ 2
      ≤ (∑ i ∈ t, (1:ℝ) ^ 2) * (∑ i ∈ t, p i ^ 2) :=
    Finset.sum_mul_sq_le_sq_mul_sq t (fun _ => (1:ℝ)) p
  simp only [one_mul, one_pow, Finset.sum_const, nsmul_eq_mul, mul_one] at hcs
  rw [← hsum]
  calc (∑ i ∈ t, p i) ^ 2
      ≤ (t.card : ℝ) * (∑ i ∈ t, p i ^ 2) := hcs
    _ ≤ (t.card : ℝ) * (∑ i ∈ s, p i ^ 2) :=
        mul_le_mul_of_nonneg_left hsq (by positivity)
