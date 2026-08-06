import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Analysis.SpecificLimits.Normed
import Mathlib.Analysis.SpecialFunctions.Exponential

/-!
# Bessel function `J₀` and its elementary positivity

mathlib has no Bessel functions.  This file introduces `besselJ0` via its power series and
proves the elementary positivity `J₀(x) > 0` for `|x| < 2` — the analytic core on which the
Riemann sign-theorem's mechanism rests (LEDGER §722): since `1/√p ≤ 1/√2 < 2` for every prime,
`J₀(1/√p) > 0`, and the diagonal onset of the θ-twisted zero sum is `μ(n) · (positive Bessel
product)`.  All results below are checked against Lean's kernel with only the trusted axioms
`{propext, Classical.choice, Quot.sound}`.
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

/-- **Bessel positivity.**  `J₀(x) > 0` for `|x| < 2`, via the alternating-series bound
`J₀(x) ≥ 1 − (x/2)²`.  In particular `J₀(1/√p) > 0` for every prime `p`. -/
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
