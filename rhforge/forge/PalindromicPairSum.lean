import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Tactic

/-!
# §1787 — Palindromic pairing sum: split a range sum into middle + symmetric pairs

Pure combinatorial rearrangement (Finset split, NOT induction — the center index shifts under induction).
For any `f : ℕ → K` and `m`,

    Σ_{i < 2m+1} f i  =  f m  +  Σ_{k < m} ( f(m+1+k) + f(m-1-k) ).

With `f i = vᵢ zⁱ` for a palindromic `v` this regroups a degree-`2m` polynomial about its middle term `zᵐ`
into symmetric pairs; combined with §1786 (`pow_pair_eq_dickson`) it gives `P(z) = zᵐ·Q(z+z⁻¹)`, feeding §1785.
EXIT0 / std-axioms.
-/

variable {M : Type*} [AddCommMonoid M]

/-- **Palindromic pairing sum.**  A `range (2m+1)` sum splits into its middle term plus symmetric pairs. -/
theorem palindromic_pair_sum (f : ℕ → M) (m : ℕ) :
    ∑ i ∈ Finset.range (2 * m + 1), f i
      = f m + ∑ k ∈ Finset.range m, (f (m + 1 + k) + f (m - 1 - k)) := by
  have hsplit : ∑ i ∈ Finset.range (2 * m + 1), f i
      = (∑ i ∈ Finset.range (m + 1), f i) + ∑ k ∈ Finset.range m, f (m + 1 + k) := by
    rw [Finset.range_eq_Ico, Finset.range_eq_Ico,
        ← Finset.sum_Ico_consecutive f (Nat.zero_le (m + 1)) (by omega : m + 1 ≤ 2 * m + 1)]
    congr 1
    rw [Finset.sum_Ico_eq_sum_range]
    apply Finset.sum_congr (by congr 1; omega) (fun k _ => rfl)
  rw [hsplit, Finset.sum_range_succ, Finset.sum_add_distrib,
      ← Finset.sum_range_reflect f m]
  ac_rfl
