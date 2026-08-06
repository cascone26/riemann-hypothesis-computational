import Mathlib.Analysis.SpecialFunctions.Trigonometric.Bounds
import Mathlib.Analysis.SpecialFunctions.Exponential
import Mathlib.Analysis.PSeries
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Order.Filter.Basic

open scoped Nat

/-- `besselJ0` via its power series (local copy for verification). -/
noncomputable def besselJ0 (x : ℝ) : ℝ := ∑' k : ℕ, (-1 : ℝ) ^ k * (x / 2) ^ (2 * k) / (k ! : ℝ) ^ 2

/-- **Sharp two-term bound on `J₀`.** For `|x| < 2`,
`|J₀(x) − (1 − (x/2)²)| ≤ (x/2)⁴/4`, i.e. `1−(x/2)²−(x/2)⁴/4 ≤ J₀(x) ≤ 1−(x/2)²+(x/2)⁴/4`.
Via the alternating-series error bound at the third term `f 2 = (x/2)⁴/4`. -/
theorem besselJ0_two_term {x : ℝ} (hx : |x| < 2) :
    |besselJ0 x - (1 - (x / 2) ^ 2)| ≤ (x / 2) ^ 4 / 4 := by
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
  have hb := alternating_series_error_bound f hfa hfs 2
  rw [← hrw] at hb
  have hs2 : (∑ i ∈ Finset.range 2, (-1 : ℝ) ^ i * f i) = 1 - y := by
    rw [Finset.sum_range_succ, Finset.sum_range_one]; simp [hf]; ring
  have hf2 : f 2 = y ^ 2 / 4 := by rw [hf]; norm_num
  rw [hs2, hf2] at hb
  have : y ^ 2 / 4 = (x / 2) ^ 4 / 4 := by rw [hy]; ring
  rw [this] at hb; exact hb

#print axioms besselJ0_two_term
