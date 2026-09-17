import IntegralRpowNeg
import Mathlib.Analysis.SumIntegralComparisons
import Mathlib.Algebra.BigOperators.Intervals

/-!
# §1289 — Integral-comparison bound `∑_{k=2}^N k^{-σ} ≤ (N^{1-σ}-1)/(1-σ)` (σ ≥ 0, σ ≠ 1)

The log-type growth bound on the truncated main sum.  Since `x ↦ x^{-σ}` is antitone on `[1,N]` (for `σ ≥ 0`),
`AntitoneOn.sum_le_integral_Ico` bounds the sum `∑_{i=1}^{N-1} (i+1)^{-σ}` (reindexed to `∑_{k=2}^N k^{-σ}`,
i.e. `Finset.Ioc 1 N`) by `∫_1^N x^{-σ} = (N^{1-σ}-1)/(1-σ)` (§1288).  Combined with the modulus bound §1287
this is the size control on ζ's finite main sum feeding §1269's conditional log-growth → sharp zero-free region.
EXIT0 / std-axioms.
-/

open Finset MeasureTheory Set intervalIntegral

theorem sum_rpow_neg_le {N : ℕ} (hN : 1 ≤ N) {σ : ℝ} (hσ0 : 0 ≤ σ) (hσ1 : σ ≠ 1) :
    ∑ k ∈ Finset.Ioc 1 N, (k : ℝ) ^ (-σ) ≤ ((N : ℝ) ^ (1 - σ) - 1) / (1 - σ) := by
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  -- x ↦ x^{-σ} is antitone on [1,N]
  have hanti : AntitoneOn (fun x : ℝ => x ^ (-σ)) (Set.Icc ((1 : ℕ) : ℝ) ((N : ℕ) : ℝ)) := by
    intro x hx y hy hxy
    refine Real.rpow_le_rpow_of_nonpos ?_ hxy (neg_nonpos.mpr hσ0)
    have : ((1 : ℕ) : ℝ) ≤ x := hx.1
    simpa using lt_of_lt_of_le one_pos (by simpa using this)
  -- sum-to-integral comparison
  have key := hanti.sum_le_integral_Ico hN
  -- reindex ∑_{i∈Ico 1 N} (i+1)^{-σ} = ∑_{k∈Ioc 1 N} k^{-σ}
  have hreindex : ∑ i ∈ Finset.Ico 1 N, ((i + 1 : ℕ) : ℝ) ^ (-σ)
      = ∑ k ∈ Finset.Ioc 1 N, (k : ℝ) ^ (-σ) := by
    rw [Finset.sum_Ico_add' (fun k : ℕ => ((k : ℝ) ^ (-σ))) 1 N 1]
    have hset : Finset.Ico (1 + 1) (N + 1) = Finset.Ioc 1 N := by
      ext x; simp only [Finset.mem_Ico, Finset.mem_Ioc]; omega
    rw [hset]
  rw [hreindex] at key
  -- the integral endpoint ↑(1:ℕ) = (1:ℝ); substitute §1288
  rw [show ((1 : ℕ) : ℝ) = (1 : ℝ) from Nat.cast_one] at key
  rw [integral_rpow_neg_Icc hNR hσ1] at key
  exact key

#print axioms sum_rpow_neg_le
