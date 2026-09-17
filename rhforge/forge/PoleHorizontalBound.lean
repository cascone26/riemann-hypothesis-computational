import Mathlib.Analysis.SpecialFunctions.Pow.Complex
import Mathlib.Analysis.SpecialFunctions.Pow.Real

set_option maxHeartbeats 1000000

open Complex

/-- §1633 — horizontal decay of the pole term: `‖(x²/2)/((σ+iT)-1)‖ ≤ (x²/2)/|T|` (`T ≠ 0`).
The pole term's contribution on a horizontal segment at height `T` decays like `1/|T|` — its half of the
rectangle's top/bottom edges vanishing. -/
theorem norm_pole_le_horizontal {x : ℝ} (hx : 0 < x) {σ T : ℝ} (hT : T ≠ 0) :
    ‖(x : ℂ) ^ (2 : ℂ) / 2 / (((σ : ℂ) + (T : ℝ) * I) - 1)‖ ≤ x ^ (2 : ℝ) / 2 / |T| := by
  have hdim : ((σ : ℂ) + (T : ℝ) * I - 1).im = T := by
    simp [Complex.add_im, Complex.mul_im, Complex.sub_im, Complex.one_im]
  have hden_ge : |T| ≤ ‖(σ : ℂ) + (T : ℝ) * I - 1‖ := by
    have := Complex.abs_im_le_norm ((σ : ℂ) + (T : ℝ) * I - 1)
    rwa [hdim] at this
  have hTpos : (0 : ℝ) < |T| := abs_pos.mpr hT
  have hnum : ‖(x : ℂ) ^ (2 : ℂ) / 2‖ = x ^ (2 : ℝ) / 2 := by
    rw [norm_div, Complex.norm_cpow_eq_rpow_re_of_pos hx]
    simp [Complex.re_ofNat]
  rw [norm_div, hnum]
  apply div_le_div_of_nonneg_left (by positivity) hTpos hden_ge
