import NormLSeriesLe
import NormZetaGeInv
import Mathlib.NumberTheory.LSeries.Dirichlet

/-!
# §1422 — closed ζ lower bound `‖ζ s‖ ≥ (∑' n, ‖term ↗μ s n‖)⁻¹` for `1 < Re s`

Chains §1421 (`‖L↗μ s‖ ≤ ∑'‖term ↗μ s n‖`) into §1420 (`‖L↗μ s‖ ≤ M ⟹ M⁻¹ ≤ ‖ζ s‖`).  Absolute summability of
the Möbius terms is free in ℂ (`summable_norm_iff` + `LSeriesSummable_moebius_iff`).  Result: a fully **explicit,
`g₀`-free** lower bound on `‖ζ s‖` — the value `∑'‖term ↗μ s n‖ = ∑ |μ n|/n^σ ≤ ζ(σ)` is a concrete constant,
supplying §1419's `‖ζ c‖`.  EXIT0 / std-axioms.
-/

open ArithmeticFunction
open scoped LSeries.notation Moebius

theorem norm_riemannZeta_ge_moebius_tsum_inv {s : ℂ} (hs : 1 < s.re) :
    (∑' n, ‖LSeries.term ↗μ s n‖)⁻¹ ≤ ‖riemannZeta s‖ := by
  have habs : Summable (fun n => ‖LSeries.term ↗μ s n‖) :=
    summable_norm_iff.mpr (LSeriesSummable_moebius_iff.mpr hs)
  exact norm_riemannZeta_ge_inv hs (norm_LSeries_le habs)

#print axioms norm_riemannZeta_ge_moebius_tsum_inv
