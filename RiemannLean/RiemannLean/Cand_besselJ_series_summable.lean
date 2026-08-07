import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Analysis.SpecificLimits.Normed
open scoped Nat
noncomputable def besselJ (n : ℕ) (x : ℝ) : ℝ :=
  ∑' k : ℕ, (-1 : ℝ) ^ k / ((k ! : ℝ) * ((k + n)! : ℝ)) * (x / 2) ^ (2 * k + n)

theorem besselJ_series_summable (n : ℕ) (x : ℝ) :
    Summable (fun k : ℕ => (-1 : ℝ) ^ k / ((k ! : ℝ) * ((k + n)! : ℝ)) * (x / 2) ^ (2 * k + n)) := by
  apply Summable.of_norm
  have hcomp := (Real.summable_pow_div_factorial ((x / 2) ^ 2)).mul_left (|x / 2| ^ n)
  refine hcomp.of_nonneg_of_le (fun k => norm_nonneg _) (fun k => ?_)
  have hkpos : (0 : ℝ) < (k ! : ℝ) := by positivity
  have hknpos : (0 : ℝ) < ((k + n)! : ℝ) := by positivity
  have hkn1 : (1 : ℝ) ≤ ((k + n)! : ℝ) := by
    exact_mod_cast Nat.one_le_iff_ne_zero.mpr (Nat.factorial_ne_zero (k + n))
  have heq : ‖(-1 : ℝ) ^ k / ((k ! : ℝ) * ((k + n)! : ℝ)) * (x / 2) ^ (2 * k + n)‖
      = |x / 2| ^ (2 * k + n) / ((k ! : ℝ) * ((k + n)! : ℝ)) := by
    rw [Real.norm_eq_abs, abs_mul, abs_div, abs_pow, abs_neg, abs_one, one_pow, abs_mul,
      abs_of_pos hkpos, abs_of_pos hknpos, abs_pow]
    ring
  rw [heq]
  calc |x / 2| ^ (2 * k + n) / ((k ! : ℝ) * ((k + n)! : ℝ))
      ≤ |x / 2| ^ (2 * k + n) / ((k ! : ℝ) * 1) := by gcongr
    _ = |x / 2| ^ n * (((x / 2) ^ 2) ^ k / (k ! : ℝ)) := by
        rw [mul_one, show 2 * k + n = n + 2 * k by ring, pow_add, pow_mul, sq_abs, mul_div_assoc]

theorem besselJ_zero_of_pos {n : ℕ} (hn : 0 < n) : besselJ n 0 = 0 := by
  unfold besselJ
  convert tsum_zero with k
  have : 2 * k + n ≠ 0 := by omega
  simp [zero_pow this]

#print axioms besselJ_series_summable
#print axioms besselJ_zero_of_pos
