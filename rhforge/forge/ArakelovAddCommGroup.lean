import Arakelov

set_option maxHeartbeats 1000000

open IsDedekindDomain NumberField

namespace ArakelovDivisor

variable (K : Type*) [Field K] [NumberField K]

/-- ℕ-scaling of an Arakelov divisor: componentwise. -/
noncomputable instance : SMul ℕ (ArakelovDivisor K) :=
  ⟨fun n D => ⟨n • D.finite, fun p => n • D.infinite p⟩⟩

/-- ℤ-scaling of an Arakelov divisor: componentwise. -/
noncomputable instance : SMul ℤ (ArakelovDivisor K) :=
  ⟨fun n D => ⟨n • D.finite, fun p => n • D.infinite p⟩⟩

/-- §1649 — `ArakelovDivisor K` is an additive commutative group. Transported along the injective
componentwise embedding into `(HeightOneSpectrum (RingOfIntegers K) →₀ ℤ) × (InfinitePlace K → ℝ)`
(both factors already `AddCommGroup`), whose zero/add/neg/sub/nsmul/zsmul all match our componentwise
instances definitionally (`rfl`). This upgrades the finite-divisor embedding `embFinite` from a mere
function into something that can be bundled as an `AddMonoidHom` (next brick). -/
noncomputable instance addCommGroup : AddCommGroup (ArakelovDivisor K) :=
  Function.Injective.addCommGroup
    (fun D => (D.finite, D.infinite))
    (fun D E h => ext K (congrArg Prod.fst h) (congrArg Prod.snd h))
    rfl (fun _ _ => rfl) (fun _ => rfl) (fun _ _ => rfl) (fun _ _ => rfl) (fun _ _ => rfl)

end ArakelovDivisor

#print axioms ArakelovDivisor.addCommGroup
