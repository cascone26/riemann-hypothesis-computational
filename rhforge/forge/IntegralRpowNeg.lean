import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

/-!
# §1288 — Real integral `∫_1^N x^{-σ} dx = (N^{1-σ}-1)/(1-σ)` (bound for the truncated main sum)

The real counterpart of §1281, for the integral-comparison bound on `∑_{k≤N} k^{-σ}` (the modulus RHS of
§1287).  With `AntitoneOn.sum_le_integral_Ico` this gives `∑_{k=2}^N k^{-σ} ≤ (N^{1-σ}-1)/(1-σ)`, which for
`σ` in the critical strip is the log-type growth of ζ's main sum → §1269.  EXIT0 / std-axioms.
-/

open MeasureTheory Set intervalIntegral

theorem integral_rpow_neg_Icc {N : ℝ} (hN : 1 ≤ N) {σ : ℝ} (hσ : σ ≠ 1) :
    ∫ x in (1 : ℝ)..N, x ^ (-σ) = (N ^ (1 - σ) - 1) / (1 - σ) := by
  have hr : (-σ : ℝ) ≠ -1 := fun h => hσ (by linarith)
  have h0 : (0 : ℝ) ∉ Set.uIcc (1 : ℝ) N := by
    rw [Set.uIcc_of_le hN, Set.mem_Icc]; push_neg; intro h; linarith
  rw [integral_rpow (Or.inr ⟨hr, h0⟩), Real.one_rpow, show (-σ + 1 : ℝ) = 1 - σ from by ring]

#print axioms integral_rpow_neg_Icc
