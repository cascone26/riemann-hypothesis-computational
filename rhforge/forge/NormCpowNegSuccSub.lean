import Mathlib.Analysis.SpecialFunctions.Pow.Deriv
import Mathlib.Analysis.Calculus.MeanValue

open Set

/-!
# §1516 — total-variation term: `‖(i+1)^{-s} - i^{-s}‖ ≤ ‖s‖·i^{-(σ+1)}` for `i ≥ 1`, `0 ≤ σ`, `s ≠ 0`

Mean-value inequality on `[i, i+1]` for `f(y) = (y:ℂ)^{-s}`, `f'(y) = -s·(y:ℂ)^{-s-1}`
(`hasDerivAt_ofReal_cpow_const`, `y ≥ 1 ≠ 0`, `-s ≠ 0`).  On the interval `‖f'(y)‖ = ‖s‖·y^{-(σ+1)} ≤ ‖s‖·i^{-(σ+1)}`
(exponent `-(σ+1) ≤ 0`, base `y ≥ i`), so `‖f(i+1)-f(i)‖ ≤ ‖s‖·i^{-(σ+1)}·1`.  Summable total-variation bound
completing the Abel estimate §1515 toward the explicit `L(χ)` growth bound.  EXIT0 / std-axioms.
-/

theorem norm_cpow_neg_succ_sub_le {s : ℂ} (hs : 0 ≤ s.re) (hs0 : s ≠ 0) {i : ℕ} (hi : 1 ≤ i) :
    ‖((i + 1 : ℕ) : ℂ) ^ (-s) - (i : ℂ) ^ (-s)‖ ≤ ‖s‖ * (i : ℝ) ^ (-(s.re + 1)) := by
  have hi0 : (0 : ℝ) < (i : ℝ) := by exact_mod_cast hi
  set f : ℝ → ℂ := fun y => (y : ℂ) ^ (-s) with hf
  set f' : ℝ → ℂ := fun y => -s * (y : ℂ) ^ (-s - 1) with hf'
  have hIcc : Convex ℝ (Icc (i : ℝ) (i + 1)) := convex_Icc _ _
  have hderiv : ∀ y ∈ Icc (i : ℝ) (i + 1), HasDerivWithinAt f (f' y) (Icc (i : ℝ) (i + 1)) y := by
    intro y hy
    have hy0 : y ≠ 0 := ne_of_gt (lt_of_lt_of_le hi0 hy.1)
    have h := hasDerivAt_ofReal_cpow_const hy0 (r := -s) (neg_ne_zero.mpr hs0)
    simpa [hf, hf', mul_comm] using h.hasDerivWithinAt
  have hbound : ∀ y ∈ Icc (i : ℝ) (i + 1), ‖f' y‖ ≤ ‖s‖ * (i : ℝ) ^ (-(s.re + 1)) := by
    intro y hy
    have hy0 : (0 : ℝ) < y := lt_of_lt_of_le hi0 hy.1
    rw [hf', norm_mul, norm_neg, Complex.norm_cpow_eq_rpow_re_of_pos hy0]
    have hre : (-s - 1).re = -(s.re + 1) := by simp; ring
    rw [hre]
    apply mul_le_mul_of_nonneg_left _ (norm_nonneg s)
    exact Real.rpow_le_rpow_of_nonpos hi0 hy.1 (by linarith)
  have hmain := hIcc.norm_image_sub_le_of_norm_hasDerivWithin_le hderiv hbound
    (left_mem_Icc.mpr (by linarith)) (right_mem_Icc.mpr (by linarith))
  have hcast : ‖f (i + 1) - f (i : ℝ)‖ = ‖((i + 1 : ℕ) : ℂ) ^ (-s) - (i : ℂ) ^ (-s)‖ := by
    rw [hf]; push_cast; ring_nf
  rw [hcast] at hmain
  calc ‖((i + 1 : ℕ) : ℂ) ^ (-s) - (i : ℂ) ^ (-s)‖
      ≤ ‖s‖ * (i : ℝ) ^ (-(s.re + 1)) * ‖(i : ℝ) + 1 - (i : ℝ)‖ := hmain
    _ = ‖s‖ * (i : ℝ) ^ (-(s.re + 1)) := by simp

#print axioms norm_cpow_neg_succ_sub_le
