import Arakelov
import ArakelovAddCommGroup

set_option maxHeartbeats 1000000

open IsDedekindDomain NumberField

namespace ArakelovDivisor

variable (K : Type*) [Field K] [NumberField K]

/-- §1650 — the finite-divisor embedding, bundled as an additive group homomorphism
`(HeightOneSpectrum (RingOfIntegers K) →₀ ℤ) →+ ArakelovDivisor K`. This packages
`embFinite_zero`/`embFinite_add` (already proven) together with the `AddCommGroup`
instance on `ArakelovDivisor K` (§1649) into the standard bundled-morphism API — the
finite Arakelov divisors now sit inside the full Arakelov group as the range of a genuine
`AddMonoidHom`, ready for `AddMonoidHom.range`/kernel-style lemmas downstream. -/
noncomputable def embFiniteHom : (HeightOneSpectrum (RingOfIntegers K) →₀ ℤ) →+ ArakelovDivisor K where
  toFun := embFinite K
  map_zero' := embFinite_zero K
  map_add' := embFinite_add K

@[simp] theorem embFiniteHom_apply (D : HeightOneSpectrum (RingOfIntegers K) →₀ ℤ) :
    embFiniteHom K D = embFinite K D := rfl

/-- The bundled embedding is injective, matching `embFinite_injective`. -/
theorem embFiniteHom_injective : Function.Injective (embFiniteHom K) :=
  embFinite_injective K

end ArakelovDivisor

#print axioms ArakelovDivisor.embFiniteHom
#print axioms ArakelovDivisor.embFiniteHom_injective
