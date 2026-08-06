import Mathlib.Analysis.SpecialFunctions.Gamma.Digamma
import Mathlib.Analysis.SpecialFunctions.Complex.Analytic

/-!
# The Weil explicit-formula quadratic form — toward the RH positivity criterion

RH is equivalent to the positivity of the Weil quadratic form
`W(φ) = archimedean(φ) + pole(φ) − prime(φ) = Σ_ρ φ̂(γ_ρ) ≥ 0` (over all suitable test functions).
This module begins a Lean formalization of the *criterion* itself (not just the consequence-side
sign theorem), built on mathlib's `Complex.digamma`.  This is the honest long-term target: a
machine-checkable RH-positivity statement whose components are verified as they are built.

Status (LEDGER §746–750): the RH-hard core is the arithmetic/prime fluctuation; the archimedean
part is the "smooth" piece (Connes' Sonin/prolate domain).  Here we start with the archimedean
weight `Ω(r) = Re[½(ψ(¼+ir/2)+ψ(¼−ir/2))] − log π`, the density against which `|φ̂|²` is integrated.
-/

open Complex

/-- The archimedean weight of the Weil explicit formula (symmetrized real form),
`Ω(r) = Re[½(ψ(¼+ir/2)+ψ(¼−ir/2))] − log π`. -/
noncomputable def archWeight (r : ℝ) : ℝ :=
  (((Complex.digamma (1 / 4 + Complex.I * r / 2)
      + Complex.digamma (1 / 4 - Complex.I * r / 2)) / 2).re) - Real.log Real.pi

/-- The archimedean weight is even in `r` (manifest from the symmetric definition). -/
theorem archWeight_even (r : ℝ) : archWeight (-r) = archWeight r := by
  unfold archWeight
  congr 2
  rw [Complex.ofReal_neg]
  ring_nf
