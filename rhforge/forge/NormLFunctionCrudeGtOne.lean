import TwistTailNormBound
import Mathlib.NumberTheory.LSeries.DirichletContinuation

open Finset Filter Topology DirichletCharacter
open scoped LSeries.notation

/-! §1522 crude LFunction bound on Re s > 1. EXIT0 / std-axioms. -/

theorem norm_LFunction_le_crude_gt_one {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M} (hχ : χ ≠ 1)
    {s : ℂ} (hs : 1 < s.re) {N : ℕ} (hN : 1 ≤ N) :
    ‖LFunction χ s‖ ≤ ((N + 1 : ℕ) : ℝ)
      + (((N + 1 : ℕ) : ℝ) ^ (-s.re) * M + ‖s‖ * (N : ℝ) ^ (-s.re) / s.re * M) := by
  have hσ : 0 < s.re := lt_trans one_pos hs
  have hmain : ‖∑ n ∈ Finset.range (N + 1), LSeries.term (↗χ) s n‖ ≤ ((N + 1 : ℕ) : ℝ) := by
    calc ‖∑ n ∈ Finset.range (N + 1), LSeries.term (↗χ) s n‖
        ≤ ∑ n ∈ Finset.range (N + 1), ‖LSeries.term (↗χ) s n‖ := norm_sum_le _ _
      _ ≤ ∑ _n ∈ Finset.range (N + 1), (1 : ℝ) := by
          refine Finset.sum_le_sum (fun n _ => ?_)
          rcases eq_or_ne n 0 with h | h
          · subst h; simp [LSeries.term]
          · rw [LSeries.term_of_ne_zero h, norm_div,
              Complex.norm_natCast_cpow_of_pos (Nat.pos_of_ne_zero h) s]
            have hn1 : (1 : ℝ) ≤ (n : ℝ) ^ s.re :=
              Real.one_le_rpow (by exact_mod_cast Nat.one_le_iff_ne_zero.mpr h) hσ.le
            rw [div_le_one (by positivity)]
            exact le_trans (χ.norm_le_one _) hn1
      _ = ((N + 1 : ℕ) : ℝ) := by simp
  have htail := twist_tail_norm_bound hχ hs hN
  rw [LFunction_eq_LSeries χ hs]
  calc ‖LSeries (↗χ) s‖
      = ‖(∑ n ∈ Finset.range (N + 1), LSeries.term (↗χ) s n)
          + (LSeries (↗χ) s - ∑ n ∈ Finset.range (N + 1), LSeries.term (↗χ) s n)‖ := by
        congr 1; ring
    _ ≤ ‖∑ n ∈ Finset.range (N + 1), LSeries.term (↗χ) s n‖
          + ‖LSeries (↗χ) s - ∑ n ∈ Finset.range (N + 1), LSeries.term (↗χ) s n‖ := norm_add_le _ _
    _ ≤ ((N + 1 : ℕ) : ℝ)
          + (((N + 1 : ℕ) : ℝ) ^ (-s.re) * M + ‖s‖ * (N : ℝ) ^ (-s.re) / s.re * M) := by
        linarith [hmain, htail]

#print axioms norm_LFunction_le_crude_gt_one
