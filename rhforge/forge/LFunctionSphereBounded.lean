import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.Analysis.Complex.CauchyIntegral

/-!
# §1478 — `LFunction χ` is bounded on every sphere (`χ ≠ 1`)

The growth-bound existence the Borel–Carathéodory factor bound §1471 needs (`hFsphere : ‖F‖ ≤ Mζ` on the sphere).
`LFunction χ` is entire (`differentiable_LFunction`), hence continuous, and the sphere `sphere c R` is compact
(`ℂ` is a `ProperSpace`), so a continuous function on it is bounded: `∃ M, ∀ s ∈ sphere c R, ‖LFunction χ s‖ ≤ M`.
Supplies the `Mζ` parameter for the L instantiation of §1471.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem exists_bound_LFunction_sphere {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    (c : ℂ) (R : ℝ) : ∃ M : ℝ, ∀ s ∈ sphere c R, ‖LFunction χ s‖ ≤ M :=
  (isCompact_sphere c R).exists_bound_of_continuousOn
    (differentiable_LFunction hχ).continuous.continuousOn

#print axioms exists_bound_LFunction_sphere
