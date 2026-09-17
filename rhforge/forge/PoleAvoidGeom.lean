import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Analysis.Normed.Field.Basic

/-!
# §1371 — pole avoidance for a disk centred at the DLVP test point `σ + iγ`

A concrete geometric input for §1355: a disk `ball (σ+iγ) R` of radius `R ≤ |γ|` never contains the pole `s = 1`,
because `dist(1, σ+iγ) = √((σ-1)² + γ²) ≥ |γ| ≥ R`.  This discharges §1355's `hball : 1 ∉ ball c R` when the
factorization/BC disk is centred at the test point (as it must be, so that §1369's centre-`c` bound lands exactly
on `logDeriv g (σ+iγ)`).  EXIT0 / std-axioms.
-/

open Complex

theorem one_notMem_ball_test {σ γ R : ℝ} (h : R ≤ |γ|) :
    (1 : ℂ) ∉ Metric.ball ((σ : ℂ) + I * (γ : ℝ)) R := by
  rw [Metric.mem_ball, Complex.dist_eq, not_lt]
  calc R ≤ |γ| := h
    _ = |((1 : ℂ) - ((σ : ℂ) + I * (γ : ℝ))).im| := by
        simp [Complex.sub_im, Complex.add_im, Complex.mul_im, Complex.I_im, Complex.I_re,
          Complex.ofReal_im, Complex.ofReal_re, abs_neg]
    _ ≤ ‖(1 : ℂ) - ((σ : ℂ) + I * (γ : ℝ))‖ := Complex.abs_im_le_norm _

#print axioms one_notMem_ball_test
