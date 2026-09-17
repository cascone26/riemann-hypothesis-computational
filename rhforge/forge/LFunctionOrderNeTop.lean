import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.Analysis.Meromorphic.Order
import Mathlib.Analysis.Complex.CauchyIntegral

/-!
# §1474 — `meromorphicOrderAt (LFunction χ) ≠ ⊤` on any connected set with a nonvanishing point

The L-analogue of the ζ brick §1330.  `MeromorphicOn.extract_zeros_poles` requires
`∀ u : U, meromorphicOrderAt (LFunction χ) u ≠ ⊤` (`LFunction χ` not locally identically zero).  On a connected
`U`, `exists_meromorphicOrderAt_ne_top_iff_forall` reduces this to a SINGLE point; at any `z₀ ∈ U` with
`LFunction χ z₀ ≠ 0`, the *entire* function `LFunction χ` (`differentiable_LFunction`, `χ ≠ 1`) is analytic and
nonvanishing, so its order is `0 ≠ ⊤`.  (Cleaner than ζ: no `z₀ ≠ 1` pole exception — `LFunction χ` is entire.)
This discharges the first character-specific hypothesis of the L divisor factorization.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Filter Topology

theorem meromorphicOrderAt_LFunction_ne_top {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) {U : Set ℂ} (hU : IsConnected U) {z₀ : ℂ} (hz₀ : z₀ ∈ U)
    (hz₀ne : LFunction χ z₀ ≠ 0) :
    ∀ u : U, meromorphicOrderAt (LFunction χ) u ≠ ⊤ := by
  have hana : AnalyticOnNhd ℂ (LFunction χ) Set.univ :=
    (differentiable_LFunction hχ).differentiableOn.analyticOnNhd isOpen_univ
  have hmero : MeromorphicOn (LFunction χ) U := hana.meromorphicOn.mono_set (Set.subset_univ U)
  refine (hmero.exists_meromorphicOrderAt_ne_top_iff_forall hU).mp ⟨⟨z₀, hz₀⟩, ?_⟩
  rw [meromorphicOrderAt_ne_top_iff_eventually_ne_zero (hmero z₀ hz₀)]
  have hcont : ContinuousAt (LFunction χ) z₀ := (differentiable_LFunction hχ z₀).continuousAt
  exact (hcont.eventually_ne hz₀ne).filter_mono nhdsWithin_le_nhds

#print axioms meromorphicOrderAt_LFunction_ne_top
