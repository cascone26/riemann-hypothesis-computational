import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.Analysis.Complex.CauchyIntegral

/-!
# §1473 — `LFunction χ` is `AnalyticOnNhd` on every ball (`χ ≠ 1`)

For a nontrivial Dirichlet character, `LFunction χ` is entire (`differentiable_LFunction`), hence differentiable
on any open ball, hence `AnalyticOnNhd` there (`DifferentiableOn.analyticOnNhd`, complex-analytic ⇔ holomorphic).
This is the holomorphy foundation the L half-ball zero-factorization needs: it feeds the generic factor-analyticity
bricks (§1383 `prod_rest_mul_analyticOnNhd`, §1346 `analyticAt_and_ne_zero_factor`) that build `g = LFunction χ / P`
for the Borel–Carathéodory bound §1471 at `F = LFunction χ`.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem LFunction_analyticOnNhd {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    (c : ℂ) (R : ℝ) : AnalyticOnNhd ℂ (LFunction χ) (ball c R) :=
  (differentiable_LFunction hχ).differentiableOn.analyticOnNhd isOpen_ball

#print axioms LFunction_analyticOnNhd
