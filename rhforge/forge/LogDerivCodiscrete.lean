import Mathlib.Analysis.Calculus.LogDeriv
import Mathlib.Analysis.Meromorphic.Divisor
import Mathlib.Analysis.Complex.Basic

/-!
# §1339 — Converter #3: codiscrete-within equality ⟹ logDeriv agree on a punctured neighborhood

The `codiscreteWithin`→local converter for the extract_zeros_poles output.  A set codiscrete within an open `U`
lies in `𝓝[≠] z` at every `z ∈ U` (`mem_codiscreteWithin_iff_forall_mem_nhdsNE`, plus `U ∈ 𝓝 z` to resolve the
`∪ Uᶜ`), so `ζ =ᶠ[codiscreteWithin U] F` gives `ζ =ᶠ[𝓝[≠] z] F`, hence (logarithmic derivative is local,
`logDeriv_congr_nhdsNE`) `logDeriv ζ =ᶠ[𝓝[≠] z] logDeriv F`.  At the DLVP test point `z = σ+iγ` (`Re > 1`, where ζ
and `F = (∏ᶠ...)•g` are both analytic and nonvanishing) this punctured equality upgrades to `logDeriv ζ z =
logDeriv F z` by continuity — the last step of the extract_zeros_poles→ζ integration.  EXIT0 / std-axioms.
-/

open Filter Topology

theorem eventuallyEq_nhdsNE_of_codiscreteWithin {ζ F : ℂ → ℂ} {U : Set ℂ} {z : ℂ}
    (hU : IsOpen U) (hz : z ∈ U) (h : ζ =ᶠ[codiscreteWithin U] F) : ζ =ᶠ[𝓝[≠] z] F := by
  filter_upwards [mem_codiscreteWithin_iff_forall_mem_nhdsNE.1 h z hz,
    nhdsWithin_le_nhds (hU.mem_nhds hz)] with x hx hxU
  exact hx.resolve_right (not_not_intro hxU)

theorem logDeriv_eventuallyEq_nhdsNE_of_codiscreteWithin {ζ F : ℂ → ℂ} {U : Set ℂ} {z : ℂ}
    (hU : IsOpen U) (hz : z ∈ U) (h : ζ =ᶠ[codiscreteWithin U] F) :
    logDeriv ζ =ᶠ[𝓝[≠] z] logDeriv F :=
  logDeriv_congr_nhdsNE (eventuallyEq_nhdsNE_of_codiscreteWithin hU hz h)

#print axioms logDeriv_eventuallyEq_nhdsNE_of_codiscreteWithin
