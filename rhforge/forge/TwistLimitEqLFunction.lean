import TwistLimit
import Mathlib.NumberTheory.LSeries.DirichletContinuation

open Finset Filter Topology DirichletCharacter
open scoped LSeries.notation

/-! §1529 twistLimit χ = LFunction χ on Re s > 1 (limit uniqueness). EXIT0 / std-axioms. -/

theorem twistLimit_eq_LFunction {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M} (hχ : χ ≠ 1)
    {s : ℂ} (hs : 1 < s.re) :
    twistLimit χ s = LFunction χ s := by
  have hσ : 0 < s.re := lt_trans one_pos hs
  have hSeq : (fun n => ∑ k ∈ Finset.Icc (1 : ℕ) n, χ (k : ZMod M) * (k : ℂ) ^ (-s))
      = (fun n => ∑ k ∈ Finset.range (n + 1), LSeries.term (↗χ) s k) := by
    funext n
    have hins : Finset.range (n + 1) = insert 0 (Finset.Icc 1 n) := by
      ext x; simp only [Finset.mem_range, Finset.mem_Icc, Finset.mem_insert]; omega
    rw [hins, Finset.sum_insert (by simp), show LSeries.term (↗χ) s 0 = 0 from by simp [LSeries.term],
      zero_add]
    apply Finset.sum_congr rfl; intro k hk
    rw [Finset.mem_Icc] at hk
    rw [LSeries.term_of_ne_zero (by omega : k ≠ 0)]
    show χ (k : ZMod M) * (k : ℂ) ^ (-s) = χ (k : ZMod M) / (k : ℂ) ^ s
    rw [Complex.cpow_neg, div_eq_mul_inv]
  have hLS : Tendsto (fun n => ∑ k ∈ Finset.Icc (1 : ℕ) n, χ (k : ZMod M) * (k : ℂ) ^ (-s)) atTop
      (𝓝 (LSeries (↗χ) s)) := by
    rw [hSeq]
    exact (LSeriesSummable_of_one_lt_re χ hs).hasSum.tendsto_sum_nat.comp (tendsto_add_atTop_nat 1)
  have huniq := tendsto_nhds_unique (tendsto_twistLimit hχ hσ) hLS
  rw [huniq, LFunction_eq_LSeries χ hs]

#print axioms twistLimit_eq_LFunction
