/-
Copyright (c) 2026 Jacob Cascone. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Jacob Cascone
-/
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Analysis.SpecificLimits.Normed
import Mathlib.Analysis.SpecialFunctions.Exponential

/-!
# Bessel functions `J₀` and `J₁` (power-series definition) and their elementary bounds

mathlib currently has no Bessel functions.  This file introduces the Bessel functions of the first
kind of orders `0` and `1` via their power series and proves elementary two-sided bounds valid below
their first positive zero:

* `besselJ0 x = ∑ₖ (-1)^k /(k!)² (x/2)^{2k}`,  `besselJ1 x = ∑ₖ (-1)^k /(k!(k+1)!) (x/2)^{2k+1}`.

## Main results

* `besselJ0_series_summable`, `besselJ0_zero` (`J₀ 0 = 1`), `besselJ0_even`.
* `besselJ0_pos`: `0 < J₀ x` for `|x| < 2` (i.e. below the first zero `j₀,₁ ≈ 2.405`).
* `besselJ0_ge_one_sub`, `besselJ0_le_one`, and the sharp two-term enclosure `besselJ0_two_term`
  (`|J₀ x − (1 − (x/2)²)| ≤ (x/2)⁴/4`) with corollary `besselJ0_le_two_term`.
* `besselJ1_pos`: `0 < J₁ x` for `0 < x`, `x² < 8` (below the first zero `j₁,₁ ≈ 3.832`).
* `besselJ1_le_half`: `J₁ x ≤ x/2` on the same range.

All bounds follow from the Leibniz alternating-series error estimate
(`alternating_series_error_bound`), the series being alternating with decreasing terms on the stated
ranges.

## TODO

Generalize to `J_ν` of arbitrary order, relate to Bessel's differential equation, and add the
integral representation.  The order-`0,1` results here are self-contained and require none of that
infrastructure.
-/

open scoped Nat

/-- The Bessel function of the first kind, order `0`, via its power series
`J₀(x) = ∑ₖ (-1)^k /(k!)² (x/2)^{2k}`. -/
noncomputable def besselJ0 (x : ℝ) : ℝ :=
  ∑' k : ℕ, (-1 : ℝ) ^ k / (k ! : ℝ) ^ 2 * (x / 2) ^ (2 * k)

/-- The `J₀` power series is (absolutely) summable, by comparison with `(x²/4)^k/k!`. -/
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

/-- `J₀(0) = 1`. -/
theorem besselJ0_zero : besselJ0 0 = 1 := by
  unfold besselJ0
  rw [tsum_eq_single 0]
  · simp
  · intro k hk
    have : (2 * k) ≠ 0 := by omega
    simp [zero_pow this]

/-- `J₀` is even: `J₀(-x) = J₀(x)` (each series term depends on `x` only through `(x/2)^{2k}`). -/
theorem besselJ0_even (x : ℝ) : besselJ0 (-x) = besselJ0 x := by
  unfold besselJ0
  congr 1; ext k
  rw [show (-x) / 2 = -(x / 2) by ring, (even_two_mul k).neg_pow]

/-- **Bessel positivity.**  `J₀(x) > 0` for `|x| < 2`, via the alternating-series bound
`J₀(x) ≥ 1 − (x/2)²`. -/
theorem besselJ0_pos {x : ℝ} (hx : |x| < 2) : 0 < besselJ0 x := by
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
    gcongr
    nlinarith [hk1]
  have hb := alternating_series_error_bound f hfa hfs 1
  rw [← hrw] at hb
  have hs1 : (∑ i ∈ Finset.range 1, (-1 : ℝ) ^ i * f i) = 1 := by
    rw [Finset.sum_range_one]; simp [hf]
  have hf1 : f 1 = y := by simp [hf]
  rw [hs1, hf1, abs_le] at hb
  nlinarith [hb.1, hy1]

/-- Quantitative Bessel lower bound: `1 − (x/2)² ≤ J₀(x)` for `|x| < 2` (alternating series). -/
theorem besselJ0_ge_one_sub {x : ℝ} (hx : |x| < 2) : 1 - (x / 2) ^ 2 ≤ besselJ0 x := by
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
  have hb := alternating_series_error_bound f hfa hfs 1
  rw [← hrw] at hb
  have hs1 : (∑ i ∈ Finset.range 1, (-1 : ℝ) ^ i * f i) = 1 := by
    rw [Finset.sum_range_one]; simp [hf]
  have hf1 : f 1 = y := by simp [hf]
  rw [hs1, hf1, abs_le] at hb
  linarith [hb.1]

/-- Upper bound `J₀(x) ≤ 1` for `|x| < 2` (alternating series ≤ first term).  With
`besselJ0_ge_one_sub`, the two-sided bound `1 − (x/2)² ≤ J₀(x) ≤ 1`. -/
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

/-- **Sharp two-term bound on `J₀`.** For `|x| < 2`, `|J₀(x) − (1 − (x/2)²)| ≤ (x/2)⁴/4`, i.e.
`1−(x/2)²−(x/2)⁴/4 ≤ J₀(x) ≤ 1−(x/2)²+(x/2)⁴/4` — a sharper two-sided enclosure than
`besselJ0_ge_one_sub`, via the alternating-series error bound at the third term `f 2 = (x/2)⁴/4`. -/
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

/-- Explicit sharp upper bound (corollary of `besselJ0_two_term`): `J₀(x) ≤ 1−(x/2)²+(x/2)⁴/4`. -/
theorem besselJ0_le_two_term {x : ℝ} (hx : |x| < 2) :
    besselJ0 x ≤ 1 - (x / 2) ^ 2 + (x / 2) ^ 4 / 4 := by
  have h := besselJ0_two_term hx
  rw [abs_le] at h; linarith [h.2]

/-- Bessel `J₁` via its power series `J₁(x)=∑ₖ(-1)^k/(k!(k+1)!)(x/2)^{2k+1}`. -/
noncomputable def besselJ1 (x : ℝ) : ℝ :=
  ∑' k : ℕ, (-1 : ℝ) ^ k / ((k ! : ℝ) * ((k + 1)! : ℝ)) * (x / 2) ^ (2 * k + 1)

/-- **Bessel `J₁` positivity**: `J₁(x) > 0` for `0 < x`, `x² < 8` (below the first positive zero of
`J₁`, `j₁,₁ ≈ 3.832`). -/
theorem besselJ1_pos {x : ℝ} (hx0 : 0 < x) (hx : x ^ 2 < 8) : 0 < besselJ1 x := by
  set y : ℝ := (x / 2) ^ 2 with hy
  have hy0 : 0 ≤ y := by rw [hy]; positivity
  have hy2 : y < 2 := by rw [hy]; nlinarith [hx]
  set f : ℕ → ℝ := fun k => y ^ k / ((k ! : ℝ) * ((k + 1)! : ℝ)) with hf
  have hfnn : ∀ k, 0 ≤ f k := fun k => by rw [hf]; positivity
  have hxhalf : (0 : ℝ) < x / 2 := by linarith
  have hfactor : besselJ1 x = (x / 2) * ∑' k, (-1 : ℝ) ^ k * f k := by
    unfold besselJ1
    rw [← tsum_mul_left]
    congr 1; ext k
    rw [hf, hy, pow_add, pow_mul, pow_one]; ring
  have hfa : Antitone f := by
    apply antitone_nat_of_succ_le; intro k
    have hfac1 : ((k + 1)! : ℝ) = ((k : ℝ) + 1) * (k ! : ℝ) := by
      rw [Nat.factorial_succ]; push_cast; ring
    have hfac2 : ((k + 2)! : ℝ) = ((k : ℝ) + 2) * ((k + 1)! : ℝ) := by
      rw [show k + 2 = (k + 1) + 1 from rfl, Nat.factorial_succ]; push_cast; ring
    have hkne : (k ! : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero k)
    have hk1ne : ((k + 1)! : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero (k + 1))
    have hstep : f (k + 1) = (y / (((k : ℝ) + 1) * ((k : ℝ) + 2))) * f k := by
      show y ^ (k + 1) / (((k + 1)! : ℝ) * ((k + 2)! : ℝ))
        = (y / (((k : ℝ) + 1) * ((k : ℝ) + 2))) * (y ^ k / ((k ! : ℝ) * ((k + 1)! : ℝ)))
      rw [hfac2, hfac1, pow_succ]; field_simp
    rw [hstep]
    have hc : y / (((k : ℝ) + 1) * ((k : ℝ) + 2)) ≤ 1 := by
      rw [div_le_one (by positivity)]; nlinarith [Nat.cast_nonneg (α := ℝ) k, hy2]
    calc (y / (((k : ℝ) + 1) * ((k : ℝ) + 2))) * f k ≤ 1 * f k :=
          mul_le_mul_of_nonneg_right hc (hfnn k)
      _ = f k := one_mul _
  have hfs : Summable f := by
    have hcomp := Real.summable_pow_div_factorial y
    apply hcomp.of_nonneg_of_le hfnn
    intro k
    have hk1 : (1 : ℝ) ≤ ((k + 1)! : ℝ) := by
      exact_mod_cast Nat.one_le_iff_ne_zero.mpr (Nat.factorial_ne_zero (k + 1))
    have hkpos : (0 : ℝ) < (k ! : ℝ) := by positivity
    show y ^ k / ((k ! : ℝ) * ((k + 1)! : ℝ)) ≤ y ^ k / (k ! : ℝ)
    gcongr; nlinarith [hk1, hkpos]
  have hb := alternating_series_error_bound f hfa hfs 1
  have hs1 : (∑ i ∈ Finset.range 1, (-1 : ℝ) ^ i * f i) = 1 := by
    rw [Finset.sum_range_one]; simp [hf]
  have hf1 : f 1 = y / 2 := by rw [hf]; norm_num [Nat.factorial]
  rw [hs1, hf1, abs_le] at hb
  have hSpos : 0 < ∑' k, (-1 : ℝ) ^ k * f k := by nlinarith [hb.1, hy2]
  rw [hfactor]; positivity

/-- Upper bound `J₁(x) ≤ x/2` for `0 < x`, `x² < 8`: since `J₁(x)=(x/2)·S` with the alternating
series `S ≤ S₀ = 1`.  With `besselJ1_pos`, completes the J₁ characterization on `(0,2√2)`. -/
theorem besselJ1_le_half {x : ℝ} (hx0 : 0 < x) (hx : x ^ 2 < 8) : besselJ1 x ≤ x / 2 := by
  set y : ℝ := (x / 2) ^ 2 with hy
  have hy0 : 0 ≤ y := by rw [hy]; positivity
  have hy2 : y < 2 := by rw [hy]; nlinarith [hx]
  set f : ℕ → ℝ := fun k => y ^ k / ((k ! : ℝ) * ((k + 1)! : ℝ)) with hf
  have hfnn : ∀ k, 0 ≤ f k := fun k => by rw [hf]; positivity
  have hxhalf : (0 : ℝ) < x / 2 := by linarith
  have hfactor : besselJ1 x = (x / 2) * ∑' k, (-1 : ℝ) ^ k * f k := by
    unfold besselJ1; rw [← tsum_mul_left]; congr 1; ext k
    rw [hf, hy, pow_add, pow_mul, pow_one]; ring
  have hfa : Antitone f := by
    apply antitone_nat_of_succ_le; intro k
    have hfac1 : ((k + 1)! : ℝ) = ((k : ℝ) + 1) * (k ! : ℝ) := by
      rw [Nat.factorial_succ]; push_cast; ring
    have hfac2 : ((k + 2)! : ℝ) = ((k : ℝ) + 2) * ((k + 1)! : ℝ) := by
      rw [show k + 2 = (k + 1) + 1 from rfl, Nat.factorial_succ]; push_cast; ring
    have hkne : (k ! : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero k)
    have hk1ne : ((k + 1)! : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero (k + 1))
    have hstep : f (k + 1) = (y / (((k : ℝ) + 1) * ((k : ℝ) + 2))) * f k := by
      show y ^ (k + 1) / (((k + 1)! : ℝ) * ((k + 2)! : ℝ))
        = (y / (((k : ℝ) + 1) * ((k : ℝ) + 2))) * (y ^ k / ((k ! : ℝ) * ((k + 1)! : ℝ)))
      rw [hfac2, hfac1, pow_succ]; field_simp
    rw [hstep]
    have hc : y / (((k : ℝ) + 1) * ((k : ℝ) + 2)) ≤ 1 := by
      rw [div_le_one (by positivity)]; nlinarith [Nat.cast_nonneg (α := ℝ) k, hy2]
    calc (y / (((k : ℝ) + 1) * ((k : ℝ) + 2))) * f k ≤ 1 * f k :=
          mul_le_mul_of_nonneg_right hc (hfnn k)
      _ = f k := one_mul _
  have hfs : Summable f := by
    have hcomp := Real.summable_pow_div_factorial y
    apply hcomp.of_nonneg_of_le hfnn
    intro k
    have hk1 : (1 : ℝ) ≤ ((k + 1)! : ℝ) := by
      exact_mod_cast Nat.one_le_iff_ne_zero.mpr (Nat.factorial_ne_zero (k + 1))
    have hkpos : (0 : ℝ) < (k ! : ℝ) := by positivity
    show y ^ k / ((k ! : ℝ) * ((k + 1)! : ℝ)) ≤ y ^ k / (k ! : ℝ)
    gcongr; nlinarith [hk1, hkpos]
  have hb := alternating_series_error_bound f hfa hfs 0
  simp only [Finset.range_zero, Finset.sum_empty, sub_zero] at hb
  have hf0 : f 0 = 1 := by simp [hf]
  rw [abs_le] at hb
  have hS : (∑' k, (-1 : ℝ) ^ k * f k) ≤ 1 := by have h2 := hb.2; rwa [hf0] at h2
  rw [hfactor]
  calc (x / 2) * ∑' k, (-1 : ℝ) ^ k * f k ≤ (x / 2) * 1 :=
        mul_le_mul_of_nonneg_left hS (le_of_lt hxhalf)
    _ = x / 2 := mul_one _
