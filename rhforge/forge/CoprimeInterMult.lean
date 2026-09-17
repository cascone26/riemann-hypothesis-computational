import IntersectionMult
import Mathlib.RingTheory.Ideal.Operations

set_option maxHeartbeats 1000000

/-- §1634 — coprime divisors don't meet: if the local ideals are coprime, `interMult I J = 0`.
Connects the naive intersection multiplicity to Mathlib's `IsCoprime` on the ideal lattice. -/
theorem interMult_eq_zero_of_isCoprime {R : Type*} [CommRing R] {I J : Ideal R}
    (h : IsCoprime I J) : interMult I J = 0 :=
  (interMult_eq_zero_iff I J).mpr (Ideal.isCoprime_iff_sup_eq.mp h)

#print axioms interMult_eq_zero_of_isCoprime
