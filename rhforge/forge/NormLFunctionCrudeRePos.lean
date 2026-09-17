import TwistLimitEqOnLFunction
import TwistTailNormBoundRePos
import TwistLimit

open Finset Filter Topology DirichletCharacter
open scoped LSeries.notation

/-! §1533 crude ‖LFunction χ‖ bound on Re s > 0 (via analytic continuation). EXIT0 / std-axioms. -/

theorem norm_LFunction_le_crude_re_pos {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M} (hχ : χ ≠ 1)
    {s : ℂ} (hσ : 0 < s.re) {N : ℕ} (hN : 1 ≤ N) :
    ‖LFunction χ s‖ ≤ (N : ℝ)
      + (((N + 1 : ℕ) : ℝ) ^ (-s.re) * M + ‖s‖ * (N : ℝ) ^ (-s.re) / s.re * M) := by
  have hEq : LFunction χ s = twistLimit χ s := (twistLimit_eqOn_LFunction hχ hσ).symm
  rw [hEq]
  have hmain : ‖∑ k ∈ Finset.Icc (1 : ℕ) N, χ (k : ZMod M) * (k : ℂ) ^ (-s)‖ ≤ (N : ℝ) := by
    calc ‖∑ k ∈ Finset.Icc (1 : ℕ) N, χ (k : ZMod M) * (k : ℂ) ^ (-s)‖
        ≤ ∑ k ∈ Finset.Icc (1 : ℕ) N, ‖χ (k : ZMod M) * (k : ℂ) ^ (-s)‖ := norm_sum_le _ _
      _ ≤ ∑ _k ∈ Finset.Icc (1 : ℕ) N, (1 : ℝ) := by
          refine Finset.sum_le_sum (fun k hk => ?_)
          rw [Finset.mem_Icc] at hk
          rw [norm_mul, Complex.norm_natCast_cpow_of_pos (by omega : 0 < k) (-s), Complex.neg_re]
          have hk1 : (1 : ℝ) ≤ (k : ℝ) ^ s.re :=
            Real.one_le_rpow (by exact_mod_cast hk.1) hσ.le
          calc ‖χ (k : ZMod M)‖ * (k : ℝ) ^ (-s.re)
              ≤ 1 * (k : ℝ) ^ (-s.re) :=
                mul_le_mul_of_nonneg_right (χ.norm_le_one _) (Real.rpow_nonneg (by positivity) _)
            _ = (k : ℝ) ^ (-s.re) := one_mul _
            _ ≤ 1 := by
                rw [Real.rpow_neg (by positivity)]
                rw [inv_le_one₀ (by positivity)]; exact hk1
      _ = (N : ℝ) := by
          rw [Finset.sum_const, Nat.card_Icc]; simp
  have htail := twist_tail_norm_bound_re_pos hχ hσ hN (tendsto_twistLimit hχ hσ)
  set P := ∑ k ∈ Finset.Icc (1 : ℕ) N, χ (k : ZMod M) * (k : ℂ) ^ (-s) with hP
  calc ‖twistLimit χ s‖
      = ‖P + (twistLimit χ s - P)‖ := by congr 1; ring
    _ ≤ ‖P‖ + ‖twistLimit χ s - P‖ := norm_add_le _ _
    _ ≤ (N : ℝ) + (((N + 1 : ℕ) : ℝ) ^ (-s.re) * M + ‖s‖ * (N : ℝ) ^ (-s.re) / s.re * M) := by
        linarith [hmain, htail]

#print axioms norm_LFunction_le_crude_re_pos
