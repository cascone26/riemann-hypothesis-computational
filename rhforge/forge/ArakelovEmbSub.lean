import Arakelov

set_option maxHeartbeats 1000000

open IsDedekindDomain NumberField

namespace ArakelovDivisor

variable (K : Type*) [Field K] [NumberField K]

/-- §1647 — the finite-divisor embedding respects subtraction: `embFinite K (D - E) = embFinite K D - embFinite K E`.
With §1645 (neg), §1644-era add/zero, `embFinite` is a full additive-group embedding of the finite divisors. -/
theorem embFinite_sub (D E : HeightOneSpectrum (RingOfIntegers K) →₀ ℤ) :
    embFinite K (D - E) = embFinite K D - embFinite K E := by
  apply ext
  · rfl
  · funext p; show (0 : ℝ) = 0 - 0; simp

#print axioms embFinite_sub
