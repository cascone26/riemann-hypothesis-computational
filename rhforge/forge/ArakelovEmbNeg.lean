import Arakelov

set_option maxHeartbeats 1000000

open IsDedekindDomain NumberField

namespace ArakelovDivisor

variable (K : Type*) [Field K] [NumberField K]

/-- §1645 — the finite-divisor embedding respects negation: `embFinite K (-D) = -(embFinite K D)`.
Completes `embFinite` as an additive group homomorphism (with `embFinite_zero`/`embFinite_add`) — the finite
Arakelov divisors sit inside the full group as a subgroup. -/
theorem embFinite_neg (D : HeightOneSpectrum (RingOfIntegers K) →₀ ℤ) :
    embFinite K (-D) = -(embFinite K D) := by
  apply ext
  · rfl
  · funext p; show (0 : ℝ) = -(0 : ℝ); simp

#print axioms embFinite_neg
