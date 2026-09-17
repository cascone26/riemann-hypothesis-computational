import NormSumIocTwistUniform
import TendstoNatRpowNegBound
import Mathlib.NumberTheory.LSeries.Dirichlet

open Finset Filter Topology
open scoped LSeries.notation

/-! §1521 twisted-tail norm bound on Re s > 1. EXIT0 / std-axioms. -/

theorem twist_tail_norm_bound {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M} (hχ : χ ≠ 1)
    {s : ℂ} (hs : 1 < s.re) {m : ℕ} (hm : 1 ≤ m) :
    ‖LSeries (↗χ) s - ∑ n ∈ Finset.range (m + 1), LSeries.term (↗χ) s n‖
      ≤ ((m + 1 : ℕ) : ℝ) ^ (-s.re) * M + ‖s‖ * (m : ℝ) ^ (-s.re) / s.re * M := by
  have hσ : 0 < s.re := lt_trans one_pos hs
  have hs0 : s ≠ 0 := by rintro rfl; rw [Complex.zero_re] at hs; linarith
  set D : ℝ := ((m + 1 : ℕ) : ℝ) ^ (-s.re) * M + ‖s‖ * (m : ℝ) ^ (-s.re) / s.re * M with hD
  have hDnn : 0 ≤ D := by
    rw [hD]; positivity
  set P := ∑ n ∈ Finset.range (m + 1), LSeries.term (↗χ) s n with hP
  have hsum : HasSum (LSeries.term (↗χ) s) (LSeries (↗χ) s) :=
    (DirichletCharacter.LSeriesSummable_of_one_lt_re χ hs).hasSum
  have hrange : Tendsto (fun k => ∑ n ∈ Finset.range (k + 1), LSeries.term (↗χ) s n)
      atTop (𝓝 (LSeries (↗χ) s)) := hsum.tendsto_sum_nat.comp (tendsto_add_atTop_nat 1)
  have hblock : ∀ k, m ≤ k → ∑ n ∈ Finset.Ioc m k, LSeries.term (↗χ) s n
      = (∑ n ∈ Finset.range (k + 1), LSeries.term (↗χ) s n) - P := by
    intro k hk
    have hIoc : Finset.Ioc m k = Finset.Ico (m + 1) (k + 1) := by
      ext x; simp only [Finset.mem_Ioc, Finset.mem_Ico]; omega
    rw [hIoc, hP, Finset.range_eq_Ico, Finset.range_eq_Ico,
      ← Finset.sum_Ico_consecutive (LSeries.term (↗χ) s) (Nat.zero_le (m + 1)) (by omega : m + 1 ≤ k + 1)]
    ring
  have htend : Tendsto (fun k => ∑ n ∈ Finset.Ioc m k, LSeries.term (↗χ) s n)
      atTop (𝓝 (LSeries (↗χ) s - P)) := by
    refine Tendsto.congr' ?_ (hrange.sub_const P)
    filter_upwards [eventually_ge_atTop m] with k hk using (hblock k hk).symm
  have hterm : ∀ k, ∑ n ∈ Finset.Ioc m k, LSeries.term (↗χ) s n
      = ∑ i ∈ Finset.Ioc m k, (i : ℂ) ^ (-s) • χ (i : ZMod M) := by
    intro k; apply Finset.sum_congr rfl; intro i hi
    rw [Finset.mem_Ioc] at hi
    rw [LSeries.term_of_ne_zero (by omega : (i : ℕ) ≠ 0)]
    show χ (i : ZMod M) / (i : ℂ) ^ s = (i : ℂ) ^ (-s) • χ (i : ZMod M)
    rw [Complex.cpow_neg, smul_eq_mul, div_eq_mul_inv]; ring
  have hnorm : Tendsto (fun k => ‖∑ n ∈ Finset.Ioc m k, LSeries.term (↗χ) s n‖)
      atTop (𝓝 ‖LSeries (↗χ) s - P‖) := htend.norm
  have hb : Tendsto (fun k : ℕ => (k : ℝ) ^ (-s.re) * M + D) atTop (𝓝 D) := by
    have := tendsto_nat_rpow_neg_bound hσ (M : ℝ) D
    simpa using this
  have hab : ∀ k, ‖∑ n ∈ Finset.Ioc m k, LSeries.term (↗χ) s n‖ ≤ (k : ℝ) ^ (-s.re) * M + D := by
    intro k
    by_cases hk : m < k
    · rw [hterm k]
      have h := norm_sum_Ioc_twist_uniform_le hχ hσ hs0 hm hk
      rw [hD]; linarith [h]
    · have hempty : Finset.Ioc m k = ∅ := by
        rw [Finset.Ioc_eq_empty]; omega
      rw [hempty, Finset.sum_empty, norm_zero]
      have : (0 : ℝ) ≤ (k : ℝ) ^ (-s.re) * M := by positivity
      linarith [hDnn]
  exact le_of_tendsto_of_tendsto' hnorm hb hab

#print axioms twist_tail_norm_bound
