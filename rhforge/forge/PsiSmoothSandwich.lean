import Mathlib.NumberTheory.Chebyshev
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

set_option maxHeartbeats 1000000

/-!
# §1604 — the de-smoothing monotone sandwich: `h·ψ(x) ≤ ∫_x^{x+h} ψ ≤ h·ψ(x+h)`

The bridge from a smoothed (`ψ₁ = ∫ ψ`) asymptotic back to `ψ` itself.  `ψ` is monotone (`psi_mono`), so on
`[x, x+h]` it is squeezed between the constants `ψ(x)` and `ψ(x+h)`; integrating (`intervalIntegral.integral_mono_on`
against `integral_const`) gives the two-sided bound.  Combined with a smoothed asymptotic `∫_x^{x+h} ψ ≈ (main term
difference)` and optimizing `h`, this yields the ψ-asymptotic (PNT with error term) from the Ψ-asymptotic.
Contour-free.  EXIT0 / std-axioms.
-/

open Chebyshev MeasureTheory

/-- Lower half of the de-smoothing sandwich: `h·ψ(x) ≤ ∫_x^{x+h} ψ`. -/
theorem psi_integral_ge {x h : ℝ} (hh : 0 < h) :
    h * ψ x ≤ ∫ t in x..(x + h), ψ t := by
  have hle : x ≤ x + h := by linarith
  have hconst : (∫ _t in x..(x + h), ψ x) = h * ψ x := by
    rw [intervalIntegral.integral_const, smul_eq_mul]; ring
  rw [← hconst]
  apply intervalIntegral.integral_mono_on hle (intervalIntegrable_const) (psi_mono.intervalIntegrable)
  intro t ht
  exact psi_mono ht.1

/-- Upper half of the de-smoothing sandwich: `∫_x^{x+h} ψ ≤ h·ψ(x+h)`. -/
theorem psi_integral_le {x h : ℝ} (hh : 0 < h) :
    (∫ t in x..(x + h), ψ t) ≤ h * ψ (x + h) := by
  have hle : x ≤ x + h := by linarith
  have hconst : (∫ _t in x..(x + h), ψ (x + h)) = h * ψ (x + h) := by
    rw [intervalIntegral.integral_const, smul_eq_mul]; ring
  rw [← hconst]
  apply intervalIntegral.integral_mono_on hle (psi_mono.intervalIntegrable) (intervalIntegrable_const)
  intro t ht
  exact psi_mono ht.2

/-- §1604 — the full de-smoothing monotone sandwich. -/
theorem psi_smooth_sandwich {x h : ℝ} (hh : 0 < h) :
    h * ψ x ≤ (∫ t in x..(x + h), ψ t) ∧ (∫ t in x..(x + h), ψ t) ≤ h * ψ (x + h) :=
  ⟨psi_integral_ge hh, psi_integral_le hh⟩

#print axioms psi_smooth_sandwich
