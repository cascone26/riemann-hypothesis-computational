import LogDerivCodiscrete
import EqOfPuncturedContinuous

/-!
# §1385 — codiscrete equality + continuity ⟹ pointwise equality

Composing §1339 (`codiscreteWithin U ⟹ 𝓝[≠] z`) with §1340 (`𝓝[≠]` equality + continuity ⟹ pointwise): two
functions agreeing on a codiscrete subset of an open `U`, both continuous at `z ∈ U`, agree at `z`.  Applied to
§1343's factorization `ζ =ᶠ[codiscreteWithin (ball c R')] (∏ᶠ)•g` at a sphere point `z ∈ sphere c R ⊆ ball c R'`
(both sides continuous there, on the pole-free ball), this upgrades the dense factorization to the pointwise
`ζ z = P_half(z)·g_full(z)` that §1378's `hfac` requires.  EXIT0 / std-axioms.
-/

open Filter Topology

theorem eq_of_codiscreteWithin_of_continuousAt {ζ F : ℂ → ℂ} {U : Set ℂ} {z : ℂ}
    (hU : IsOpen U) (hz : z ∈ U) (h : ζ =ᶠ[codiscreteWithin U] F)
    (hζ : ContinuousAt ζ z) (hF : ContinuousAt F z) : ζ z = F z :=
  eq_of_eventuallyEq_nhdsNE_of_continuousAt
    (eventuallyEq_nhdsNE_of_codiscreteWithin hU hz h) hζ hF

#print axioms eq_of_codiscreteWithin_of_continuousAt
