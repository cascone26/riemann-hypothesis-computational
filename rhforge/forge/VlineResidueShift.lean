import RectangleResidue

set_option maxHeartbeats 1000000

open Complex intervalIntegral

/-!
# §1660 — THE VERTICAL-LINE RESIDUE SHIFT: crossing a simple pole deposits `2πc`

Rearrangement of §1658 into contour-shift form: for a simple pole `c·(z-p)⁻¹` and two vertical
truncated lines `Re = σ₁ < p.re` and `Re = σ₀ > p.re` (heights `[a,b]` straddling `p.im`), the
difference of the two line integrals equals the residue deposit `2πc` plus the two horizontal
connector integrals. With `c = x²/2`, `p = 1`: THE Perron pole-term shift across `s = 1` — as
`T → ∞` the connectors vanish (§1635/§1638 shape) and the deposit `2π·x²/2` is the PNT main term
(the `1/(2πi)`-normalized Perron integral picks up exactly `x²/2·(1/i)·...` bookkeeping downstream).
-/

theorem vline_diff_residue (c p : ℂ) {σ₁ σ₀ a b : ℝ}
    (h1 : σ₁ < p.re) (h2 : p.re < σ₀) (h3 : a < p.im) (h4 : p.im < b) :
    (∫ y in a..b, c * ((σ₀ : ℂ) + y * Complex.I - p)⁻¹) -
      (∫ y in a..b, c * ((σ₁ : ℂ) + y * Complex.I - p)⁻¹)
      = 2 * (Real.pi : ℂ) * c +
        Complex.I * ((∫ t in σ₁..σ₀, c * ((t : ℂ) + a * Complex.I - p)⁻¹) -
          (∫ t in σ₁..σ₀, c * ((t : ℂ) + b * Complex.I - p)⁻¹)) := by
  have hz_re : ((σ₁ : ℂ) + (a : ℝ) * Complex.I).re = σ₁ := by simp
  have hz_im : ((σ₁ : ℂ) + (a : ℝ) * Complex.I).im = a := by simp
  have hw_re : ((σ₀ : ℂ) + (b : ℝ) * Complex.I).re = σ₀ := by simp
  have hw_im : ((σ₀ : ℂ) + (b : ℝ) * Complex.I).im = b := by simp
  have h := integral_boundary_rect_residue c p ((σ₁ : ℂ) + a * Complex.I)
    ((σ₀ : ℂ) + b * Complex.I)
    (by rw [hz_re]; exact h1) (by rw [hw_re]; exact h2)
    (by rw [hz_im]; exact h3) (by rw [hw_im]; exact h4)
  rw [hz_re, hz_im, hw_re, hw_im] at h
  simp only [smul_eq_mul] at h
  linear_combination (-Complex.I) * h +
    ((∫ y in a..b, c * ((σ₀ : ℂ) + y * Complex.I - p)⁻¹) -
      (∫ y in a..b, c * ((σ₁ : ℂ) + y * Complex.I - p)⁻¹) -
      2 * (Real.pi : ℂ) * c) * Complex.I_sq

#print axioms vline_diff_residue
