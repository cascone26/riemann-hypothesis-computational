import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Analysis.SpecificLimits.Normed
import Mathlib.Analysis.SpecialFunctions.Exponential

open scoped Nat

/-- Bessel function of the first kind, integer order `n`, via its power series
`J_n(x) = ∑ₖ (-1)^k / (k!(k+n)!) (x/2)^{2k+n}`. -/
noncomputable def besselJ (n : ℕ) (x : ℝ) : ℝ :=
  ∑' k : ℕ, (-1 : ℝ) ^ k / ((k ! : ℝ) * ((k + n)! : ℝ)) * (x / 2) ^ (2 * k + n)

/-- **General Bessel positivity.**  `J_n(x) > 0` for `0 < x` and `x² < 4(n+1)` (below the first
positive zero of `J_n`).  Unifies and generalizes `J₀ > 0` on `|x|<2` (`n=0`) and `J₁ > 0` on
`x²<8` (`n=1`), via the alternating-series bound `J_n(x) ≥ (x/2)ⁿ/n! · (1 − (x/2)²/(n+1))`. -/
theorem besselJ_pos {n : ℕ} {x : ℝ} (hx0 : 0 < x) (hx : x ^ 2 < 4 * (n + 1)) :
    0 < besselJ n x := by
  set y : ℝ := (x / 2) ^ 2 with hy
  have hy0 : 0 ≤ y := by rw [hy]; positivity
  have hyn : y < (n : ℝ) + 1 := by rw [hy]; nlinarith [hx]
  set f : ℕ → ℝ := fun k => y ^ k / ((k ! : ℝ) * ((k + n)! : ℝ)) with hf
  have hfnn : ∀ k, 0 ≤ f k := fun k => by rw [hf]; positivity
  have hxhalf : (0 : ℝ) < (x / 2) ^ n := by positivity
  have hfactor : besselJ n x = (x / 2) ^ n * ∑' k, (-1 : ℝ) ^ k * f k := by
    unfold besselJ
    rw [← tsum_mul_left]
    congr 1; ext k
    rw [hf, hy, pow_add, pow_mul]; ring
  have hfa : Antitone f := by
    apply antitone_nat_of_succ_le; intro k
    have hfac1 : ((k + 1)! : ℝ) = ((k : ℝ) + 1) * (k ! : ℝ) := by
      rw [Nat.factorial_succ]; push_cast; ring
    have hfac2 : ((k + 1 + n)! : ℝ) = ((k : ℝ) + 1 + n) * ((k + n)! : ℝ) := by
      rw [show k + 1 + n = (k + n) + 1 by ring, Nat.factorial_succ]; push_cast; ring
    have hkne : (k ! : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero k)
    have hknne : ((k + n)! : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero (k + n))
    have hstep : f (k + 1) = (y / (((k : ℝ) + 1) * ((k : ℝ) + 1 + n))) * f k := by
      show y ^ (k + 1) / (((k + 1)! : ℝ) * ((k + 1 + n)! : ℝ))
        = (y / (((k : ℝ) + 1) * ((k : ℝ) + 1 + n))) * (y ^ k / ((k ! : ℝ) * ((k + n)! : ℝ)))
      rw [hfac1, hfac2, pow_succ]; field_simp
    rw [hstep]
    have hc : y / (((k : ℝ) + 1) * ((k : ℝ) + 1 + n)) ≤ 1 := by
      rw [div_le_one (by positivity)]
      nlinarith [Nat.cast_nonneg (α := ℝ) k, Nat.cast_nonneg (α := ℝ) n, hyn]
    calc (y / (((k : ℝ) + 1) * ((k : ℝ) + 1 + n))) * f k ≤ 1 * f k :=
          mul_le_mul_of_nonneg_right hc (hfnn k)
      _ = f k := one_mul _
  have hfs : Summable f := by
    have hcomp := Real.summable_pow_div_factorial y
    apply hcomp.of_nonneg_of_le hfnn
    intro k
    have hkn1 : (1 : ℝ) ≤ ((k + n)! : ℝ) := by
      exact_mod_cast Nat.one_le_iff_ne_zero.mpr (Nat.factorial_ne_zero (k + n))
    have hkpos : (0 : ℝ) < (k ! : ℝ) := by positivity
    show y ^ k / ((k ! : ℝ) * ((k + n)! : ℝ)) ≤ y ^ k / (k ! : ℝ)
    gcongr; nlinarith [hkn1, hkpos]
  have hb := alternating_series_error_bound f hfa hfs 1
  have hs1 : (∑ i ∈ Finset.range 1, (-1 : ℝ) ^ i * f i) = 1 / (n ! : ℝ) := by
    rw [Finset.sum_range_one]; simp [hf]
  rw [hs1] at hb
  have hnfac : (0 : ℝ) < (n ! : ℝ) := by positivity
  have hn1fac : (0 : ℝ) < ((n + 1)! : ℝ) := by positivity
  have hSpos : 0 < ∑' k, (-1 : ℝ) ^ k * f k := by
    rw [abs_le] at hb
    have hlow := hb.1
    have : (1 : ℝ) / (n ! : ℝ) - f 1 ≤ ∑' k, (-1 : ℝ) ^ k * f k := by linarith [hlow]
    have hf1val : f 1 = y / ((n + 1)! : ℝ) := by
      show y ^ 1 / ((1 ! : ℝ) * ((1 + n)! : ℝ)) = y / ((n + 1)! : ℝ)
      rw [Nat.factorial_one, show 1 + n = n + 1 from by ring, pow_one]
      norm_num
    have hfacrel : ((n + 1)! : ℝ) = ((n : ℝ) + 1) * (n ! : ℝ) := by
      rw [Nat.factorial_succ]; push_cast; ring
    have hpos : 0 < (1 : ℝ) / (n ! : ℝ) - f 1 := by
      rw [hf1val, hfacrel]
      rw [div_sub_div _ _ (ne_of_gt hnfac) (by positivity), div_pos_iff]
      left
      constructor
      · nlinarith [hyn, hnfac]
      · positivity
    linarith
  rw [hfactor]
  positivity

#print axioms besselJ_pos
