import Psi1CDeriv
import PsiSmoothSandwich
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic

set_option maxHeartbeats 1000000

/-!
# §1605 — the ℂ↔ℝ bridge: `Ψℂ(x) = ↑(Ψℝ(x))` and the smoothed-difference identity

`psi1C` (the ℂ-valued smoothed Chebyshev function fed to the Mellin inversion) is just the `Complex.ofReal` of the
real `psi1R x = ∫₁ˣ ψ` (`intervalIntegral.integral_ofReal`).  Its increments are genuine `ψ`-integrals
(`integral_add_adjacent_intervals`, `ψ` monotone ⟹ interval-integrable), so the de-smoothing sandwich §1604 applies
verbatim to the increments of the Perron output §1601/§1602:  `h·ψ(x) ≤ Ψℝ(x+h) - Ψℝ(x) ≤ h·ψ(x+h)`.  This is the
final gluing lemma converting the contour-shift asymptotic of `Ψ` into the `ψ`-asymptotic (PNT with error term).
EXIT0 / std-axioms.
-/

open Chebyshev MeasureTheory

/-- The real smoothed Chebyshev function `Ψℝ(x) = ∫₁ˣ ψ`. -/
noncomputable def psi1R (x : ℝ) : ℝ := ∫ t in (1 : ℝ)..x, ψ t

/-- Bridge: the ℂ-valued `psi1C` is the `ofReal` of the real `psi1R`. -/
theorem psi1C_eq_ofReal (x : ℝ) : psi1C x = ((psi1R x : ℝ) : ℂ) := by
  rw [psi1C, psi1R]
  exact intervalIntegral.integral_ofReal

/-- Increment identity: `Ψℝ(x+h) - Ψℝ(x) = ∫_x^{x+h} ψ`. -/
theorem psi1R_diff (x h : ℝ) :
    psi1R (x + h) - psi1R x = ∫ t in x..(x + h), ψ t := by
  rw [psi1R, psi1R]
  rw [← intervalIntegral.integral_add_adjacent_intervals
    (a := (1:ℝ)) (b := x) (c := x + h) psi_mono.intervalIntegrable psi_mono.intervalIntegrable]
  ring

/-- §1605 — the de-smoothing sandwich, applied to the smoothed Chebyshev increments. -/
theorem psi1R_sandwich {x h : ℝ} (hh : 0 < h) :
    h * ψ x ≤ psi1R (x + h) - psi1R x ∧ psi1R (x + h) - psi1R x ≤ h * ψ (x + h) := by
  rw [psi1R_diff]
  exact psi_smooth_sandwich hh

#print axioms psi1R_sandwich
