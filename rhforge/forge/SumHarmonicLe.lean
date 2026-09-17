import Mathlib.Analysis.SumIntegralComparisons
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Algebra.BigOperators.Intervals

/-!
# §1307 — Harmonic bound `∑_{k=2}^N k⁻¹ ≤ log N` (the `σ = 1` case of §1289)

At `σ = 1` the integral-comparison bound §1289 breaks (`integral_rpow` excludes `r = -1`; the divided form
`(N^{1-σ}-1)/(1-σ)` is `0/0`).  The correct bound is the harmonic one: `x ↦ x^{-1}` antitone on `[1,N]`,
`AntitoneOn.sum_le_integral_Ico`, and `∫_1^N x⁻¹ = log N` (`integral_inv_of_pos`).  This covers the `σ = 1`
vertical line that a Borel–Carathéodory disk centred at `1 + it₀` passes through.  EXIT0 / std-axioms.
-/

open Finset MeasureTheory Set intervalIntegral

theorem sum_rpow_neg_one_le {N : ℕ} (hN : 1 ≤ N) :
    ∑ k ∈ Finset.Ioc 1 N, (k : ℝ) ^ (-1 : ℝ) ≤ Real.log N := by
  have hanti : AntitoneOn (fun x : ℝ => x ^ (-1 : ℝ)) (Set.Icc ((1 : ℕ) : ℝ) ((N : ℕ) : ℝ)) := by
    intro x hx y hy hxy
    refine Real.rpow_le_rpow_of_nonpos ?_ hxy (by norm_num)
    have : ((1 : ℕ) : ℝ) ≤ x := hx.1
    simpa using lt_of_lt_of_le one_pos (by simpa using this)
  have key := hanti.sum_le_integral_Ico hN
  have hreindex : ∑ i ∈ Finset.Ico 1 N, ((i + 1 : ℕ) : ℝ) ^ (-1 : ℝ)
      = ∑ k ∈ Finset.Ioc 1 N, (k : ℝ) ^ (-1 : ℝ) := by
    rw [Finset.sum_Ico_add' (fun k : ℕ => ((k : ℝ) ^ (-1 : ℝ))) 1 N 1]
    have hset : Finset.Ico (1 + 1) (N + 1) = Finset.Ioc 1 N := by
      ext x; simp only [Finset.mem_Ico, Finset.mem_Ioc]; omega
    rw [hset]
  rw [hreindex] at key
  rw [show ((1 : ℕ) : ℝ) = (1 : ℝ) from Nat.cast_one] at key
  have hint : ∫ x in (1 : ℝ)..(N : ℝ), x ^ (-1 : ℝ) = Real.log N := by
    have hcongr : ∀ x ∈ Set.uIcc (1 : ℝ) (N : ℝ), x ^ (-1 : ℝ) = x⁻¹ :=
      fun x _ => Real.rpow_neg_one x
    rw [intervalIntegral.integral_congr hcongr,
      integral_inv_of_pos one_pos (by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hN), div_one]
  rw [hint] at key
  exact key

#print axioms sum_rpow_neg_one_le
