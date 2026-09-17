import SmoothedPerronDisplay
import Mathlib.MeasureTheory.Group.Integral

set_option maxHeartbeats 2000000

open Complex Chebyshev MeasureTheory

/-!
# §1699 — THE Ψ-PERRON BRIDGE: `2π·Ψ(x) = ∫_ℝ F(σ₀+it) dt` in the working `F`-spelling

Connects the banked Mellin-inversion display (§1602: `Ψ(x) = (1/2π)∫ x^{σ₀+1-iy}·(-ζ'/ζ)/(s(s+1))`
at `s = σ₀-iy`) to the exact integrand spelling used throughout the winding/Perron arcs
(`F(σ₀+it) = -ζ'/ζ(σ₀+it)·x^{s+1}/(s(s+1))`): the display integrand IS `F(σ₀-iy)`, so the
substitution `y ↦ -t` (volume is negation-invariant) gives `2π·Ψ(x) = ∫_ℝ F(σ₀+it) dt` — the
full-line integral whose truncations §1698 controls. The last conventions-bridge of the chain.
-/

theorem two_pi_psi1C_eq_perron_integral {σ₀ x : ℝ} (hσ₀ : 1 < σ₀) (hx : 0 < x) :
    (2 * Real.pi : ℝ) • psi1C x = ∫ t : ℝ,
      -logDeriv riemannZeta ((σ₀ : ℝ) + (t : ℝ) * Complex.I) *
          (x : ℂ) ^ (((σ₀ : ℝ) + (t : ℝ) * Complex.I) + 1) /
        (((σ₀ : ℝ) + (t : ℝ) * Complex.I) * (((σ₀ : ℝ) + (t : ℝ) * Complex.I) + 1)) := by
  have hπ : (2 * Real.pi : ℝ) ≠ 0 := by
    have := Real.pi_pos
    positivity
  rw [smoothed_perron_display hσ₀ hx, smul_smul, mul_one_div, div_self hπ, one_smul]
  -- pointwise: the display integrand at y equals F(σ₀ + i(-y))
  have hpt : ∀ y : ℝ,
      (x : ℂ) ^ (-(((-σ₀ - 1 : ℝ) : ℂ) + (y : ℝ) * I)) •
        (-logDeriv riemannZeta ((σ₀ : ℂ) - (y : ℝ) * I)
          / (((σ₀ : ℂ) - (y : ℝ) * I) * (((σ₀ : ℂ) - (y : ℝ) * I) + 1)))
      = (fun t : ℝ =>
          -logDeriv riemannZeta ((σ₀ : ℝ) + (t : ℝ) * Complex.I) *
              (x : ℂ) ^ (((σ₀ : ℝ) + (t : ℝ) * Complex.I) + 1) /
            (((σ₀ : ℝ) + (t : ℝ) * Complex.I) *
              (((σ₀ : ℝ) + (t : ℝ) * Complex.I) + 1))) (-y) := by
    intro y
    simp only []
    have hs : ((σ₀ : ℝ) : ℂ) + ((-y : ℝ) : ℂ) * Complex.I = (σ₀ : ℂ) - (y : ℝ) * I := by
      push_cast
      ring
    have hexp : -(((-σ₀ - 1 : ℝ) : ℂ) + (y : ℝ) * I) = ((σ₀ : ℂ) - (y : ℝ) * I) + 1 := by
      push_cast
      ring
    rw [hs, hexp, smul_eq_mul]
    ring
  calc (∫ y : ℝ, (x : ℂ) ^ (-(((-σ₀ - 1 : ℝ) : ℂ) + (y : ℝ) * I)) •
        (-logDeriv riemannZeta ((σ₀ : ℂ) - (y : ℝ) * I)
          / (((σ₀ : ℂ) - (y : ℝ) * I) * (((σ₀ : ℂ) - (y : ℝ) * I) + 1))))
      = ∫ y : ℝ, (fun t : ℝ =>
          -logDeriv riemannZeta ((σ₀ : ℝ) + (t : ℝ) * Complex.I) *
              (x : ℂ) ^ (((σ₀ : ℝ) + (t : ℝ) * Complex.I) + 1) /
            (((σ₀ : ℝ) + (t : ℝ) * Complex.I) *
              (((σ₀ : ℝ) + (t : ℝ) * Complex.I) + 1))) (-y) := by
        congr 1
        funext y
        exact hpt y
    _ = ∫ t : ℝ,
          -logDeriv riemannZeta ((σ₀ : ℝ) + (t : ℝ) * Complex.I) *
              (x : ℂ) ^ (((σ₀ : ℝ) + (t : ℝ) * Complex.I) + 1) /
            (((σ₀ : ℝ) + (t : ℝ) * Complex.I) *
              (((σ₀ : ℝ) + (t : ℝ) * Complex.I) + 1)) := by
        have h := integral_neg_eq_self (fun t : ℝ =>
          -logDeriv riemannZeta ((σ₀ : ℝ) + (t : ℝ) * Complex.I) *
              (x : ℂ) ^ (((σ₀ : ℝ) + (t : ℝ) * Complex.I) + 1) /
            (((σ₀ : ℝ) + (t : ℝ) * Complex.I) *
              (((σ₀ : ℝ) + (t : ℝ) * Complex.I) + 1))) (volume : Measure ℝ)
        convert h using 2

#print axioms two_pi_psi1C_eq_perron_integral
