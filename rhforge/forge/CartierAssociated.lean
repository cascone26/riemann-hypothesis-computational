import IntersectionMult
import Mathlib.RingTheory.Ideal.Span

set_option maxHeartbeats 1000000

/-- §1643 — the local Cartier intersection number is invariant under associates: if `f` and `f'` differ by a unit
(`Associated f f'`), then `(f)·(g) = (f')·(g)`.  Well-definedness on principal divisors of an integral ring
(they depend only on the ideal `span {f}`, not the generator). -/
theorem cartierInterMult_of_associated {R : Type*} [CommRing R] [IsDomain R] {f f' g : R}
    (h : Associated f f') : cartierInterMult f g = cartierInterMult f' g := by
  unfold cartierInterMult
  rw [Ideal.span_singleton_eq_span_singleton.mpr h]

#print axioms cartierInterMult_of_associated
