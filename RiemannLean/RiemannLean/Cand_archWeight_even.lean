import Mathlib.Analysis.SpecialFunctions.Gamma.Digamma
import Mathlib.Analysis.SpecialFunctions.Complex.Analytic
open Complex
/-- The archimedean weight of the Weil explicit formula (symmetrized real form),
`Ω(r) = Re[½(ψ(¼+ir/2)+ψ(¼−ir/2))] − log π`.  Building block toward formalizing the Weil
positivity criterion (mathlib has `Complex.digamma`). -/
noncomputable def archWeight (r : ℝ) : ℝ :=
  (((Complex.digamma (1/4 + Complex.I * r / 2)
      + Complex.digamma (1/4 - Complex.I * r / 2)) / 2).re) - Real.log Real.pi
/-- The archimedean weight is even in `r` (manifest from the symmetric definition). -/
theorem archWeight_even (r : ℝ) : archWeight (-r) = archWeight r := by
  unfold archWeight
  congr 2
  rw [Complex.ofReal_neg]
  ring_nf

#print axioms archWeight_even
