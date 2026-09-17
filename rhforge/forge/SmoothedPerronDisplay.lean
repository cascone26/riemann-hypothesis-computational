import SmoothedPerron
import Psi1CMellinForm
import Mathlib.Analysis.MellinTransform
import Mathlib.Analysis.MellinInversion

set_option maxHeartbeats 1000000

/-!
# §1602 — the SMOOTHED PERRON DISPLAY corollary: the explicit vertical integral

Unfolding `mellinInv` in the smoothed Perron formula §1601 and rewriting `mellin Ψ` on the vertical line via §1598
(`mellin_psi1C_line_eq`) gives the fully explicit representation of the smoothed Chebyshev function as a contour
integral of `-ζ'/ζ`:

  `Ψ(x) = (1/2π) • ∫_ℝ x^{-((-σ₀-1)+iy)} • (-ζ'/ζ(σ₀-iy))/((σ₀-iy)(σ₀-iy+1)) dy`.

This is the display form fed to the contour-shift main-term extraction (residue at `s=1` = `x²/2`) and the error
control from the zero-free region §1575.  EXIT0 / std-axioms.
-/

open Complex Chebyshev MeasureTheory

/-- The Mellin transform of `Ψ` on the line `Re = -σ₀-1`, in explicit `-ζ'/ζ` form (per-point §1598). -/
theorem mellin_psi1C_line_eq {σ₀ : ℝ} (hσ₀ : 1 < σ₀) (y : ℝ) :
    mellin psi1C (((-σ₀ - 1 : ℝ) : ℂ) + (y : ℝ) * I)
      = -logDeriv riemannZeta ((σ₀ : ℂ) - (y : ℝ) * I)
          / (((σ₀ : ℂ) - (y : ℝ) * I) * (((σ₀ : ℂ) - (y : ℝ) * I) + 1)) := by
  have hsyre : ((σ₀ : ℂ) - (y : ℝ) * I).re = σ₀ := by
    simp [Complex.sub_re, Complex.mul_re]
  have hsyne : (σ₀ : ℂ) - (y : ℝ) * I ≠ 0 := by
    intro h; have := congrArg Complex.re h; rw [hsyre] at this; simp at this; linarith
  have hsy1ne : (σ₀ : ℂ) - (y : ℝ) * I + 1 ≠ 0 := by
    intro h; have := congrArg Complex.re h
    rw [Complex.add_re, hsyre, Complex.one_re] at this; simp at this; linarith
  have hsygt : 1 < ((σ₀ : ℂ) - (y : ℝ) * I).re := by rw [hsyre]; exact hσ₀
  have harg : ((-σ₀ - 1 : ℝ) : ℂ) + (y : ℝ) * I = -((σ₀ : ℂ) - (y : ℝ) * I) - 1 := by
    push_cast; ring
  rw [harg]
  have h98 := mellin_psi1C_eq hsygt
  rw [eq_div_iff (mul_ne_zero hsyne hsy1ne)]
  linear_combination h98

/-- §1602 — the explicit vertical-integral display of the smoothed Perron formula. -/
theorem smoothed_perron_display {σ₀ x : ℝ} (hσ₀ : 1 < σ₀) (hx : 0 < x) :
    psi1C x = (1 / (2 * Real.pi)) • ∫ y : ℝ,
      (x : ℂ) ^ (-(((-σ₀ - 1 : ℝ) : ℂ) + (y : ℝ) * I)) •
        (-logDeriv riemannZeta ((σ₀ : ℂ) - (y : ℝ) * I)
          / (((σ₀ : ℂ) - (y : ℝ) * I) * (((σ₀ : ℂ) - (y : ℝ) * I) + 1))) := by
  have hfun : (fun y : ℝ => (x : ℂ) ^ (-(((-σ₀ - 1 : ℝ) : ℂ) + (y : ℝ) * I)) •
        mellin psi1C (((-σ₀ - 1 : ℝ) : ℂ) + (y : ℝ) * I))
      = (fun y : ℝ => (x : ℂ) ^ (-(((-σ₀ - 1 : ℝ) : ℂ) + (y : ℝ) * I)) •
        (-logDeriv riemannZeta ((σ₀ : ℂ) - (y : ℝ) * I)
          / (((σ₀ : ℂ) - (y : ℝ) * I) * (((σ₀ : ℂ) - (y : ℝ) * I) + 1)))) := by
    funext y; rw [mellin_psi1C_line_eq hσ₀ y]
  rw [← smoothed_perron hσ₀ hx]
  simp only [mellinInv]
  rw [hfun]

#print axioms smoothed_perron_display
