/-
  Basic verification that Mathlib RH infrastructure works.
-/

import Mathlib.NumberTheory.LSeries.RiemannZeta
import Mathlib.NumberTheory.ArithmeticFunction
import Mathlib.NumberTheory.Harmonic.EulerMascheroni

-- Verify key definitions exist
#check RiemannHypothesis
#check riemannZeta
#check ArithmeticFunction.sigma
#check Real.eulerMascheroniConstant
