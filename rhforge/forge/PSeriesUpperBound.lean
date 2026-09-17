import Mathlib.Analysis.SumIntegralComparisons
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Analysis.PSeries

set_option maxHeartbeats 1000000

open Finset intervalIntegral

/-!
# §1685 — the p-series EXPLICIT BOUND: `Σ' 1/n^σ ≤ 1 + 1/(σ-1)` for `σ > 1`

The quantitative upper companion of `Real.summable_one_div_nat_rpow`, by the integral test
(`AntitoneOn.sum_le_integral` + `integral_rpow`). At `σ = 1 + 1/log T` this is `≤ 1 + log T` —
turning §1678's `‖ζ(c)‖ ≥ (Σ 1/n^σ)⁻¹` into the concrete `‖ζ(c)‖ ≥ 1/(1 + log T)`, the last
numeric input of the box-edge `O(log²T)`.
-/

theorem tsum_one_div_nat_rpow_le {σ : ℝ} (hσ : 1 < σ) :
    ∑' n : ℕ, 1 / (n : ℝ) ^ σ ≤ 1 + 1 / (σ - 1) := by
  have hσ0 : (0 : ℝ) < σ := by linarith
  apply Real.tsum_le_of_sum_range_le (fun n => by positivity)
  intro N₀
  calc ∑ n ∈ Finset.range N₀, 1 / (n : ℝ) ^ σ
      ≤ ∑ n ∈ Finset.range (N₀ + 2), 1 / (n : ℝ) ^ σ :=
        Finset.sum_le_sum_of_subset_of_nonneg
          (by intro x hx; rw [Finset.mem_range] at hx ⊢; omega)
          (fun n _ _ => by positivity)
    _ ≤ 1 + 1 / (σ - 1) := ?_
  set N : ℕ := N₀ with hNdef
  clear_value N
  -- split off n = 0, 1
  rw [Finset.sum_range_succ', Finset.sum_range_succ']
  have h0 : 1 / ((0 : ℕ) : ℝ) ^ σ = 0 := by
    simp [Real.zero_rpow hσ0.ne']
  have h1 : 1 / (((0 : ℕ) + 1 : ℕ) : ℝ) ^ σ = 1 := by
    norm_num
  rw [h0, h1]
  -- the tail via the integral test
  have htail : ∑ i ∈ Finset.range N, 1 / (((i + 1 : ℕ) + 1 : ℕ) : ℝ) ^ σ
      ≤ ∫ x in (1 : ℝ)..(1 + N), x ^ (-σ) := by
    have hanti : AntitoneOn (fun x : ℝ => x ^ (-σ)) (Set.Icc 1 (1 + (N : ℝ))) := by
      intro x hx y hy hxy
      have hx1 : (1 : ℝ) ≤ x := hx.1
      have hxpos : (0 : ℝ) < x := by linarith
      have hypos : (0 : ℝ) < y := by linarith [hy.1]
      simp only [Real.rpow_neg hxpos.le, Real.rpow_neg hypos.le]
      rw [inv_eq_one_div, inv_eq_one_div]
      exact one_div_le_one_div_of_le (Real.rpow_pos_of_pos hxpos σ)
        (Real.rpow_le_rpow hxpos.le hxy hσ0.le)
    have h := hanti.sum_le_integral (a := N)
    calc ∑ i ∈ Finset.range N, 1 / (((i + 1 : ℕ) + 1 : ℕ) : ℝ) ^ σ
        = ∑ i ∈ Finset.range N, ((1 : ℝ) + ((i + 1 : ℕ) : ℝ)) ^ (-σ) := by
          apply Finset.sum_congr rfl
          intro i _
          rw [Real.rpow_neg (x := (1 : ℝ) + ((i + 1 : ℕ) : ℝ)) (by positivity),
            inv_eq_one_div]
          congr 2
          push_cast
          ring
      _ ≤ ∫ x in (1 : ℝ)..(1 + N), x ^ (-σ) := h
  -- evaluate the integral
  have hint : (∫ x in (1 : ℝ)..(1 + N), x ^ (-σ)) ≤ 1 / (σ - 1) := by
    have hne : (-σ) ≠ -1 := by
      intro h
      apply hσ.ne'
      linarith [neg_injective h]
    have hmem : (0 : ℝ) ∉ Set.uIcc (1 : ℝ) (1 + N) := by
      rw [Set.uIcc_of_le (by have := Nat.cast_nonneg (α := ℝ) N; linarith)]
      intro h
      have := (Set.mem_Icc.mp h).1
      linarith
    rw [integral_rpow (Or.inr ⟨hne, hmem⟩)]
    have h1p : (1 : ℝ) ^ (-σ + 1) = 1 := Real.one_rpow _
    rw [h1p, show (-σ + 1) = -(σ - 1) by ring, div_neg, ← neg_div, neg_sub]
    have hX : (0 : ℝ) ≤ (1 + (N : ℝ)) ^ (-(σ - 1)) :=
      Real.rpow_nonneg (by positivity) _
    gcongr <;> linarith
  linarith [htail, hint]

#print axioms tsum_one_div_nat_rpow_le
