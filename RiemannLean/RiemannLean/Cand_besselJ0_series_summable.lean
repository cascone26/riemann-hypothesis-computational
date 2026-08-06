import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Analysis.SpecialFunctions.Exponential
open scoped Nat
/-- Summability of the Bessel J₀ series `∑ (-1)^k /(k!)² (x/2)^{2k}` (first step toward
    Bessel J₀ in mathlib, which lacks Bessel functions entirely). By comparison with the
    exponential series `(x²/4)^k/k!` using `(k!)² ≥ k!`. -/
theorem besselJ0_series_summable (x : ℝ) :
    Summable (fun k : ℕ => (-1 : ℝ) ^ k / (k ! : ℝ) ^ 2 * (x / 2) ^ (2 * k)) := by
  apply Summable.of_norm
  have hcomp := Real.summable_pow_div_factorial ((x / 2) ^ 2)
  refine hcomp.of_nonneg_of_le (fun k => norm_nonneg _) (fun k => ?_)
  have hkpos : (0 : ℝ) < (k ! : ℝ) := by positivity
  have hk1 : (1 : ℝ) ≤ (k ! : ℝ) := by
    exact_mod_cast Nat.one_le_iff_ne_zero.mpr (Nat.factorial_ne_zero k)
  rw [pow_mul]
  have hnorm : ‖(-1 : ℝ) ^ k / (k ! : ℝ) ^ 2 * ((x / 2) ^ 2) ^ k‖
      = ((x / 2) ^ 2) ^ k / (k ! : ℝ) ^ 2 := by
    rw [norm_mul, norm_div, norm_pow, norm_neg, norm_one, one_pow, norm_pow, norm_pow]
    rw [Real.norm_eq_abs, abs_of_pos hkpos, Real.norm_eq_abs,
        abs_of_nonneg (by positivity : (0 : ℝ) ≤ (x / 2) ^ 2)]
    ring
  rw [hnorm]
  gcongr
  nlinarith [hk1]

#print axioms besselJ0_series_summable
