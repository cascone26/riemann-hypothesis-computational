import CoprimeInterMult
import Mathlib.RingTheory.Ideal.Operations

set_option maxHeartbeats 1000000

/-- §1636 — coprime elements give disjoint principal divisors: if `f, g` are coprime in `R`
then their local Cartier intersection number vanishes, `(f)·(g) = 0`.  The element-level form of §1634,
via `isCoprime_span_singleton_iff`. -/
theorem cartierInterMult_eq_zero_of_isCoprime {R : Type*} [CommRing R] {f g : R}
    (h : IsCoprime f g) : cartierInterMult f g = 0 :=
  interMult_eq_zero_of_isCoprime ((Ideal.isCoprime_span_singleton_iff f g).mpr h)

#print axioms cartierInterMult_eq_zero_of_isCoprime
