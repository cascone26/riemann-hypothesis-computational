import Mathlib.Analysis.Complex.Basic
import Mathlib.Topology.Separation.Basic

/-!
# §1340 — Punctured-neighborhood equality + continuity ⟹ equality at the point

If `f =ᶠ[𝓝[≠] z] g` and both are continuous at `z` (in ℂ, where `𝓝[≠] z` is `NeBot`), then `f z = g z`: both
functions tend to their values at `z` along `𝓝[≠] z`, and along that filter `g` inherits `f`'s limit, so the two
limits coincide (`tendsto_nhds_unique`).  This upgrades §1339's punctured `logDeriv` equality to the pointwise value
`logDeriv ζ z = logDeriv F z` at the DLVP test point `z = σ+iγ` (`Re > 1`), where both log-derivatives are continuous
(ζ and `F = (∏ᶠ...)•g` analytic and nonvanishing there).  EXIT0 / std-axioms.
-/

open Filter Topology

theorem eq_of_eventuallyEq_nhdsNE_of_continuousAt {f g : ℂ → ℂ} {z : ℂ}
    (hfg : f =ᶠ[𝓝[≠] z] g) (hf : ContinuousAt f z) (hg : ContinuousAt g z) : f z = g z := by
  have h1 : Tendsto g (𝓝[≠] z) (𝓝 (f z)) :=
    (hf.tendsto.mono_left nhdsWithin_le_nhds).congr' hfg
  have h2 : Tendsto g (𝓝[≠] z) (𝓝 (g z)) :=
    hg.tendsto.mono_left nhdsWithin_le_nhds
  exact tendsto_nhds_unique h1 h2

#print axioms eq_of_eventuallyEq_nhdsNE_of_continuousAt
