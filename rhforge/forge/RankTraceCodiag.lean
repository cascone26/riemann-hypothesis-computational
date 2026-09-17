import Mathlib

/-!
# Co-diagonal rank–trace inequality (RH session 2026-08-29)

The full Alpöge–Furman rank–trace inequality (AF eq (L)) in the simultaneously-diagonalizable case:
for a PSD spectrum `p ≥ 0` and any spectrum `q`,
`rank P = #{pᵢ≠0} ≥ 2 tr P + 4 tr Q − 4·npl(Q) − ‖P+Q‖²_HS = 2∑p + 4∑q − 4#{q>0} − ∑(pᵢ+qᵢ)²`.
Per eigenvalue (using `pᵢ≥0`) the integrality shadows `(p-1)²≥0`, `(q-2)²≥0` absorb the cross term.
Combines `RankTraceSingle` (P-part) and `RankTraceQ` (Q-part). Verified: 0 violations / 200k random.
-/

open Finset

/-- Co-diagonal rank–trace: `2∑p + 4∑q − 4#{q>0} − ∑(pᵢ+qᵢ)² ≤ #{pᵢ≠0}` for `p ≥ 0`. -/
theorem rank_trace_codiagonal {ι : Type*} (s : Finset ι) (p q : ι → ℝ)
    (hp : ∀ i ∈ s, 0 ≤ p i) :
    2 * (∑ i ∈ s, p i) + 4 * (∑ i ∈ s, q i)
        - 4 * (s.filter (fun i => 0 < q i)).card - (∑ i ∈ s, (p i + q i) ^ 2)
      ≤ (s.filter (fun i => p i ≠ 0)).card := by
  classical
  have per : ∀ i ∈ s,
      2 * p i + 4 * q i - 4 * (if 0 < q i then (1:ℝ) else 0) - (p i + q i) ^ 2
        ≤ (if p i ≠ 0 then (1:ℝ) else 0) := by
    intro i hi
    have hpi := hp i hi
    by_cases hq : 0 < q i
    · rw [if_pos hq]
      by_cases hpz : p i = 0
      · rw [if_neg (not_not.mpr hpz)]
        nlinarith [sq_nonneg (q i - 2), hpi]
      · rw [if_pos hpz]
        nlinarith [sq_nonneg (p i - 1), sq_nonneg (q i - 2), hpi, hq.le]
    · rw [if_neg hq]
      have hqle : q i ≤ 0 := not_lt.mp hq
      by_cases hpz : p i = 0
      · rw [if_neg (not_not.mpr hpz)]
        nlinarith [sq_nonneg (p i + q i), hpi, hqle, sq_nonneg (q i)]
      · rw [if_pos hpz]
        nlinarith [sq_nonneg (p i - 1), sq_nonneg (q i), hpi, hqle]
  have hsum : ∑ i ∈ s,
      (2 * p i + 4 * q i - 4 * (if 0 < q i then (1:ℝ) else 0) - (p i + q i) ^ 2)
        ≤ ∑ i ∈ s, (if p i ≠ 0 then (1:ℝ) else 0) := Finset.sum_le_sum per
  have hRq : ∑ i ∈ s, (if 0 < q i then (1:ℝ) else 0)
      = ((s.filter (fun i => 0 < q i)).card : ℝ) := by
    rw [Finset.sum_ite, Finset.sum_const, Finset.sum_const_zero, add_zero,
        nsmul_eq_mul, mul_one]
  have hRp : ∑ i ∈ s, (if p i ≠ 0 then (1:ℝ) else 0)
      = ((s.filter (fun i => p i ≠ 0)).card : ℝ) := by
    rw [Finset.sum_ite, Finset.sum_const, Finset.sum_const_zero, add_zero,
        nsmul_eq_mul, mul_one]
  have hL : ∑ i ∈ s,
      (2 * p i + 4 * q i - 4 * (if 0 < q i then (1:ℝ) else 0) - (p i + q i) ^ 2)
      = 2 * (∑ i ∈ s, p i) + 4 * (∑ i ∈ s, q i)
          - 4 * (∑ i ∈ s, (if 0 < q i then (1:ℝ) else 0)) - (∑ i ∈ s, (p i + q i) ^ 2) := by
    simp only [Finset.sum_sub_distrib, Finset.sum_add_distrib, Finset.mul_sum]
  rw [hL, hRq, hRp] at hsum
  exact hsum
