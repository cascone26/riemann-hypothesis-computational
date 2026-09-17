import Mathlib.Analysis.SpecialFunctions.Pow.Asymptotics
import Mathlib.Order.Filter.AtTopBot.Archimedean

open Filter Topology

/-!
# §1520 — vanishing of the `n^{-σ}` block term: `Tendsto (fun n => n^{-σ}·M + D) atTop (𝓝 D)`

For `σ > 0`, `(n:ℝ)^{-σ} → 0` (`tendsto_rpow_neg_atTop` composed with `Nat.cast → ∞`), so the uniform block bound
§1519 (whose only `n`-dependence is the `n^{-σ}·M` term) converges, as `n → ∞`, to the tail bound `D`.  This is the
sequence whose limit `le_of_tendsto` will compare against the twisted-tail norm — the bridge from the finite block
bound §1519 to the twisted-tail estimate.  EXIT0 / std-axioms.
-/

theorem tendsto_nat_rpow_neg_bound {σ : ℝ} (hσ : 0 < σ) (M D : ℝ) :
    Tendsto (fun n : ℕ => (n : ℝ) ^ (-σ) * M + D) atTop (𝓝 (0 * M + D)) := by
  have h1 : Tendsto (fun n : ℕ => (n : ℝ) ^ (-σ)) atTop (𝓝 0) :=
    (tendsto_rpow_neg_atTop hσ).comp tendsto_natCast_atTop_atTop
  exact ((h1.mul_const M).add_const D)

#print axioms tendsto_nat_rpow_neg_bound
