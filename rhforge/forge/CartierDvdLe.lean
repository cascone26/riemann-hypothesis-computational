import IntersectionMult
import Mathlib.RingTheory.Ideal.Span

set_option maxHeartbeats 1000000

/-- §1641 — the local Cartier intersection number is MONOTONE under divisibility: if `f ∣ f'`
then `(f)·(g) ≤ (f')·(g)`.  A larger divisor (more vanishing) meets `g` at least as much.  Via
`span {f'} ≤ span {f}` (from `f ∣ f'`) and antitonicity of `interMult` in the first ideal (§1629). -/
theorem cartierInterMult_dvd_le {R : Type*} [CommRing R] {f f' g : R} (h : f ∣ f') :
    cartierInterMult f g ≤ cartierInterMult f' g :=
  interMult_antitone_left (Ideal.span {g}) (Ideal.span_singleton_le_span_singleton.mpr h)

#print axioms cartierInterMult_dvd_le
