import IntersectionMult
import Mathlib.RingTheory.Ideal.Operations

set_option maxHeartbeats 1000000

/-- §1639 — vanishing Cartier intersection number CHARACTERIZES coprimality:
`(f)·(g) = 0 ↔ IsCoprime f g`.  The intersection number is zero exactly when the two principal divisors share
no common zero — the geometric content, now an iff. -/
theorem cartierInterMult_eq_zero_iff_isCoprime {R : Type*} [CommRing R] (f g : R) :
    cartierInterMult f g = 0 ↔ IsCoprime f g := by
  rw [cartierInterMult_eq_zero_iff, ← Ideal.isCoprime_iff_sup_eq, Ideal.isCoprime_span_singleton_iff]

#print axioms cartierInterMult_eq_zero_iff_isCoprime
