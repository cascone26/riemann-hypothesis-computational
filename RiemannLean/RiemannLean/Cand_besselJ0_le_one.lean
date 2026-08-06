import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Analysis.SpecificLimits.Normed
import Mathlib.Analysis.SpecialFunctions.Exponential
open scoped Nat
noncomputable def besselJ0 (x : ℝ) : ℝ := ∑' k : ℕ, (-1 : ℝ) ^ k / (k ! : ℝ) ^ 2 * (x / 2) ^ (2 * k)
/-- Upper bound `J₀(x) ≤ 1` for `|x| < 2` (alternating series ≤ first term). With
`besselJ0_ge_one_sub` gives the two-sided `1−(x/2)² ≤ J₀(x) ≤ 1`. -/
theorem besselJ0_le_one {x : ℝ} (hx : |x| < 2) : besselJ0 x ≤ 1 := by
  set y : ℝ := (x / 2) ^ 2 with hy
  have hy0 : 0 ≤ y := by rw [hy]; positivity
  have hy1 : y < 1 := by rw [hy]; rw [abs_lt] at hx; nlinarith [hx.1, hx.2]
  set f : ℕ → ℝ := fun k => y ^ k / (k ! : ℝ) ^ 2 with hf
  have hfnn : ∀ k, 0 ≤ f k := fun k => by rw [hf]; positivity
  have hrw : besselJ0 x = ∑' k, (-1 : ℝ) ^ k * f k := by
    unfold besselJ0; congr 1; ext k; rw [hf, hy, pow_mul]; ring
  have hfa : Antitone f := by
    apply antitone_nat_of_succ_le; intro k
    have hfac : ((k + 1)! : ℝ) = ((k : ℝ) + 1) * (k ! : ℝ) := by
      rw [Nat.factorial_succ]; push_cast; ring
    have hkne : (k ! : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero k)
    have hk1ne : ((k : ℝ) + 1) ≠ 0 := by positivity
    have hstep : f (k + 1) = (y / ((k : ℝ) + 1) ^ 2) * f k := by
      show y ^ (k + 1) / ((k + 1)! : ℝ) ^ 2 = (y / ((k : ℝ) + 1) ^ 2) * (y ^ k / (k ! : ℝ) ^ 2)
      rw [hfac, pow_succ]; field_simp
    rw [hstep]
    have hc : y / ((k : ℝ) + 1) ^ 2 ≤ 1 := by
      rw [div_le_one (by positivity)]; nlinarith [Nat.cast_nonneg (α := ℝ) k, hy1]
    calc (y / ((k : ℝ) + 1) ^ 2) * f k ≤ 1 * f k :=
          mul_le_mul_of_nonneg_right hc (hfnn k)
      _ = f k := one_mul _
  have hfs : Summable f := by
    have hcomp := Real.summable_pow_div_factorial y
    apply hcomp.of_nonneg_of_le hfnn
    intro k
    have hk1 : (1 : ℝ) ≤ (k ! : ℝ) := by
      exact_mod_cast Nat.one_le_iff_ne_zero.mpr (Nat.factorial_ne_zero k)
    show y ^ k / (k ! : ℝ) ^ 2 ≤ y ^ k / (k ! : ℝ)
    gcongr; nlinarith [hk1]
  have hb := alternating_series_error_bound f hfa hfs 0
  rw [← hrw] at hb
  simp only [Finset.range_zero, Finset.sum_empty, sub_zero] at hb
  have hf0 : f 0 = 1 := by simp [hf]
  rw [abs_le] at hb
  linarith [hb.2, hf0.symm.le, hf0.ge]

#print axioms besselJ0_le_one
