import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Algebra.Order.BigOperators.Group.Finset

/-!
# §1438 — crude uniform partial-sum bound `∑_{k∈Ioc 1 N} k^{-σ} ≤ N` for `0 ≤ σ`

Each term `k^{-σ} ≤ 1` (base `k ≥ 1`, exponent `-σ ≤ 0`) and there are `N-1 ≤ N` terms.  Unlike the sharp
`(N^{1-σ}-1)/(1-σ)` bound (§`sum_rpow_neg_le`), this is UNIFORM across `σ = 1` — exactly what a sphere crossing
`Re s = 1` needs for an explicit ζ growth bound.  EXIT0 / std-axioms.
-/

open Finset

theorem sum_rpow_neg_le_card {N : ℕ} {σ : ℝ} (hσ : 0 ≤ σ) :
    ∑ k ∈ Ioc 1 N, (k : ℝ) ^ (-σ) ≤ (N : ℝ) := by
  calc ∑ k ∈ Ioc 1 N, (k : ℝ) ^ (-σ)
      ≤ ∑ _k ∈ Ioc 1 N, (1 : ℝ) := by
        apply Finset.sum_le_sum
        intro k hk
        have hk1 : (1 : ℝ) ≤ (k : ℝ) := by
          have := (Finset.mem_Ioc.mp hk).1
          exact_mod_cast Nat.one_le_of_lt this
        exact Real.rpow_le_one_of_one_le_of_nonpos hk1 (by linarith)
    _ = ((Ioc 1 N).card : ℝ) := by rw [Finset.sum_const, nsmul_eq_mul, mul_one]
    _ ≤ (N : ℝ) := by rw [Nat.card_Ioc]; exact_mod_cast Nat.sub_le N 1

#print axioms sum_rpow_neg_le_card
