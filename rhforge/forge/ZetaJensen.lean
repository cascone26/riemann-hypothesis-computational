import Mathlib.Analysis.Complex.JensenFormula
import Mathlib.NumberTheory.LSeries.RiemannZeta
import Mathlib.Analysis.Meromorphic.Divisor
import ExplicitFormula

/-!
# §1361 — Jensen's formula for the Riemann zeta function (zero-counting foundation)

Specialising Mathlib's `MeromorphicOn.circleAverage_log_norm` (Jensen's formula) to `ζ`, which is meromorphic
on every closed ball.  The circle average of `log ‖ζ‖` decomposes into the divisor-weighted zero/pole counting
term `∑ᶠ u, divisor ζ u · log(R/‖c-u‖)` plus the trailing-coefficient value.  This is the foundation of the
factor bound `‖ζ/P‖ ≤ M'`: since every zero `ρ` in the open ball contributes `divisor ρ · log(R/‖c-ρ‖) ≥ 0`,
the zero count is controlled by the (boundable) circle average of `log ‖ζ‖`, which discharges §1355's `hgB`
via §1352.  EXIT0 / std-axioms.
-/

open Complex Metric Real

theorem riemannZeta_circleAverage_log_norm {c : ℂ} {R : ℝ} (hR : R ≠ 0) :
    circleAverage (fun s => Real.log ‖riemannZeta s‖) c R
      = ∑ᶠ u, (MeromorphicOn.divisor riemannZeta (closedBall c |R|) u : ℝ)
            * Real.log (R * ‖c - u‖⁻¹)
        + (MeromorphicOn.divisor riemannZeta (closedBall c |R|) c : ℝ) * Real.log R
        + Real.log ‖meromorphicTrailingCoeffAt riemannZeta c‖ := by
  have hmero : MeromorphicOn riemannZeta (closedBall c |R|) :=
    meromorphicOn_riemannZeta.mono_set (Set.subset_univ _)
  exact hmero.circleAverage_log_norm hR

#print axioms riemannZeta_circleAverage_log_norm
