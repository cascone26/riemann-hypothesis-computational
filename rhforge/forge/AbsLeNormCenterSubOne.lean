import Mathlib.Analysis.SpecialFunctions.Complex.Circle

/-!
# §1446 — `|γ| ≤ ‖(σ+iγ) - 1‖`

The imaginary part of `(σ+iγ) - 1` is `γ`, and `|Im z| ≤ ‖z‖`.  This discharges §1440's pole-distance condition
`R < ‖c-1‖` (for `c = σK+iγ`) from the scalar `R < |γ|`, and matches §1371's pole-avoidance geometry.
EXIT0 / std-axioms.
-/

open Complex

theorem abs_le_norm_center_sub_one {σ γ : ℝ} :
    |γ| ≤ ‖((σ : ℝ) : ℂ) + I * (γ : ℝ) - 1‖ := by
  have him : (((σ : ℝ) : ℂ) + I * (γ : ℝ) - 1).im = γ := by
    simp [Complex.add_im, Complex.sub_im, Complex.mul_im, Complex.one_im]
  calc |γ| = |(((σ : ℝ) : ℂ) + I * (γ : ℝ) - 1).im| := by rw [him]
    _ ≤ ‖((σ : ℝ) : ℂ) + I * (γ : ℝ) - 1‖ := Complex.abs_im_le_norm _

#print axioms abs_le_norm_center_sub_one
