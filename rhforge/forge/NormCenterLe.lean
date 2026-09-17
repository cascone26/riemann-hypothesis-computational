import Mathlib.Analysis.SpecialFunctions.Complex.Circle

/-!
# §1450 — `‖(σ+iγ)‖ ≤ |σ| + |γ|`

Triangle inequality: `‖σ + iγ‖ ≤ ‖σ‖ + ‖iγ‖ = |σ| + |γ|`.  Bounds the `‖c‖` appearing in the BC scalar of §1449
by an explicit `|σK| + |γ| ≍ |γ|`, a step toward the constant choice.  EXIT0 / std-axioms.
-/

open Complex

theorem norm_center_le {σ γ : ℝ} : ‖((σ : ℝ) : ℂ) + I * (γ : ℝ)‖ ≤ |σ| + |γ| := by
  calc ‖((σ : ℝ) : ℂ) + I * (γ : ℝ)‖
      ≤ ‖((σ : ℝ) : ℂ)‖ + ‖I * (γ : ℝ)‖ := norm_add_le _ _
    _ = |σ| + |γ| := by
        rw [Complex.norm_real, Real.norm_eq_abs, norm_mul, Complex.norm_I, one_mul,
          Complex.norm_real, Real.norm_eq_abs]

#print axioms norm_center_le
