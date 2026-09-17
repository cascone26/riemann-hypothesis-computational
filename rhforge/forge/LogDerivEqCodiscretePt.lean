import LogDerivCodiscrete
import EqOfPuncturedContinuous

/-!
# §1341 — Pointwise `logDeriv ζ z = logDeriv F z` from codiscrete equality + continuity

Combining §1339 (codiscrete-within equality ⟹ `logDeriv ζ =ᶠ[𝓝[≠] z] logDeriv F`) with §1340 (punctured equality
+ continuity ⟹ equality at the point).  For `ζ =ᶠ[codiscreteWithin U] F` on an open `U ∋ z`, if both
log-derivatives are continuous at `z`, then `logDeriv ζ z = logDeriv F z`.  This is the complete converter #3 for
the extract_zeros_poles output: at the DLVP test point `z = σ+iγ` (`Re > 1`, ζ and `F = (∏ᶠ...)•g` analytic and
nonvanishing, so both `logDeriv` continuous) it gives the pointwise value the decomposition needs.  EXIT0 / std-axioms.
-/

open Filter Topology

theorem logDeriv_eq_of_codiscreteWithin_of_continuousAt {ζ F : ℂ → ℂ} {U : Set ℂ} {z : ℂ}
    (hU : IsOpen U) (hz : z ∈ U) (h : ζ =ᶠ[codiscreteWithin U] F)
    (hζc : ContinuousAt (logDeriv ζ) z) (hFc : ContinuousAt (logDeriv F) z) :
    logDeriv ζ z = logDeriv F z :=
  eq_of_eventuallyEq_nhdsNE_of_continuousAt
    (logDeriv_eventuallyEq_nhdsNE_of_codiscreteWithin hU hz h) hζc hFc

#print axioms logDeriv_eq_of_codiscreteWithin_of_continuousAt
