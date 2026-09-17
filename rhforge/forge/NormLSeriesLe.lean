import Mathlib.NumberTheory.LSeries.Basic

/-!
# §1421 — `‖L f s‖ ≤ ∑' n, ‖term f s n‖` (norm of an L-series ≤ sum of term norms)

`LSeries f s = ∑' n, term f s n` by definition, so the triangle inequality for `tsum`
(`norm_tsum_le_tsum_norm`) gives the bound directly, given absolute summability of the terms.  Instantiated at
`f = ↗μ` (with the comparison `|μ n|/n^σ ≤ 1/n^σ`) this yields the explicit `M` for §1420's ζ lower bound.
EXIT0 / std-axioms.
-/

theorem norm_LSeries_le {f : ℕ → ℂ} {s : ℂ}
    (hsum : Summable (fun n => ‖LSeries.term f s n‖)) :
    ‖LSeries f s‖ ≤ ∑' n, ‖LSeries.term f s n‖ :=
  norm_tsum_le_tsum_norm hsum

#print axioms norm_LSeries_le
