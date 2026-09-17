import Mathlib.NumberTheory.LSeries.RiemannZeta
import Mathlib.Analysis.Meromorphic.Order
import ExplicitFormula

/-!
# §1330 — `meromorphicOrderAt ζ ≠ ⊤` on any connected set with a nonvanishing point

The Borel–Carathéodory factorization `extract_zeros_poles` requires `∀ u : U, meromorphicOrderAt ζ u ≠ ⊤`
(ζ not locally identically zero).  On a connected `U`, `exists_meromorphicOrderAt_ne_top_iff_forall` reduces this
to a SINGLE point — and at any `z₀ ∈ U` with `z₀ ≠ 1` and `ζ(z₀) ≠ 0`, ζ is analytic and nonvanishing, so its
order is `0 ≠ ⊤`.  This discharges one of the two ζ-specific hypotheses of the divisor factorization.
EXIT0 / std-axioms.
-/

open Complex Filter Topology

theorem meromorphicOrderAt_riemannZeta_ne_top {U : Set ℂ} (hU : IsConnected U)
    {z₀ : ℂ} (hz₀ : z₀ ∈ U) (hz₀1 : z₀ ≠ 1) (hz₀ne : riemannZeta z₀ ≠ 0) :
    ∀ u : U, meromorphicOrderAt riemannZeta u ≠ ⊤ := by
  have hmero : MeromorphicOn riemannZeta U := meromorphicOn_riemannZeta.mono_set (Set.subset_univ U)
  refine (hmero.exists_meromorphicOrderAt_ne_top_iff_forall hU).mp ⟨⟨z₀, hz₀⟩, ?_⟩
  rw [meromorphicOrderAt_ne_top_iff_eventually_ne_zero (hmero z₀ hz₀)]
  have hcont : ContinuousAt riemannZeta z₀ := (differentiableAt_riemannZeta hz₀1).continuousAt
  exact (hcont.eventually_ne hz₀ne).filter_mono nhdsWithin_le_nhds

#print axioms meromorphicOrderAt_riemannZeta_ne_top
