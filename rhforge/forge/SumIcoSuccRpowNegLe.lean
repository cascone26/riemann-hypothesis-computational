import Mathlib.Analysis.SumIntegralComparisons
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

open Finset

/-!
# §1517 — summed total-variation bound: `∑_{i∈Ico m N} (i+1)^{-(σ+1)} ≤ m^{-σ}/σ` (`σ>0`, `m≥1`)

Integral comparison (`AntitoneOn.sum_le_integral_Ico`) for the antitone `t ↦ t^{-(σ+1)}`:
`∑_{i∈Ico m N} (i+1)^{-(σ+1)} ≤ ∫_m^N t^{-(σ+1)} dt = (m^{-σ}-N^{-σ})/σ ≤ m^{-σ}/σ`.  This convergent tail bound,
fed with §1516's difference bound, controls the total variation in the Abel estimate §1515.  EXIT0 / std-axioms.
-/

theorem sum_Ico_succ_rpow_neg_le {σ : ℝ} (hσ : 0 < σ) {m N : ℕ} (hm : 1 ≤ m) (hmN : m ≤ N) :
    ∑ i ∈ Finset.Ico m N, ((i : ℝ) + 1) ^ (-(σ + 1)) ≤ (m : ℝ) ^ (-σ) / σ := by
  have hm0 : (0 : ℝ) < (m : ℝ) := by exact_mod_cast hm
  have hmNr : (m : ℝ) ≤ (N : ℝ) := by exact_mod_cast hmN
  have hanti : AntitoneOn (fun t : ℝ => t ^ (-(σ + 1))) (Set.Icc (m : ℝ) N) := by
    intro x hx y hy hxy
    exact Real.rpow_le_rpow_of_nonpos (lt_of_lt_of_le hm0 hx.1) hxy (by linarith)
  have hsum := hanti.sum_le_integral_Ico hmN
  have hlhs : ∑ i ∈ Finset.Ico m N, (fun t : ℝ => t ^ (-(σ + 1))) ↑(i + 1)
      = ∑ i ∈ Finset.Ico m N, ((i : ℝ) + 1) ^ (-(σ + 1)) := by
    apply Finset.sum_congr rfl; intro i _; push_cast; ring_nf
  rw [hlhs] at hsum
  have hint : ∫ x in (m : ℝ)..(N : ℝ), x ^ (-(σ + 1)) = ((N : ℝ) ^ (-σ) - (m : ℝ) ^ (-σ)) / (-σ) := by
    rw [integral_rpow (Or.inr ⟨by linarith, by
      simp only [Set.uIcc_of_le hmNr, Set.mem_Icc, not_and, not_le]
      intro h; linarith⟩)]
    norm_num
  rw [hint] at hsum
  have hNσ : (0 : ℝ) ≤ (N : ℝ) ^ (-σ) := Real.rpow_nonneg (by positivity) _
  have heq : ((N : ℝ) ^ (-σ) - (m : ℝ) ^ (-σ)) / (-σ) = ((m : ℝ) ^ (-σ) - (N : ℝ) ^ (-σ)) / σ := by
    rw [div_neg, ← neg_div, neg_sub]
  rw [heq] at hsum
  have hbound : ((m : ℝ) ^ (-σ) - (N : ℝ) ^ (-σ)) / σ ≤ (m : ℝ) ^ (-σ) / σ := by
    gcongr
    linarith [hNσ]
  linarith

#print axioms sum_Ico_succ_rpow_neg_le
