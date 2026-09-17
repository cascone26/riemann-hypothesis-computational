import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.Analysis.Meromorphic.Divisor
import Mathlib.Analysis.Complex.CauchyIntegral

/-!
# §1475 — `divisor (LFunction χ) (ball c R)` has finite support (`χ ≠ 1`)

The second character-specific hypothesis of `MeromorphicOn.extract_zeros_poles`, L-analogue of ζ §1331.  The
entire function `LFunction χ` (`differentiable_LFunction`) is meromorphic on the compact `closedBall c R`
(`ℂ` is a `ProperSpace`), so `divisor_ball_support_finite` gives finitely many zeros on the open ball — the
divisor has finite support.  With §1474 (`order ≠ ⊤`) this discharges BOTH inputs of the L divisor
factorization `LFunction χ = P·g`.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem divisor_LFunction_ball_support_finite {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) (c : ℂ) (R : ℝ) :
    (MeromorphicOn.divisor (LFunction χ) (ball c R)).support.Finite := by
  have hana : AnalyticOnNhd ℂ (LFunction χ) Set.univ :=
    (differentiable_LFunction hχ).differentiableOn.analyticOnNhd isOpen_univ
  exact (hana.meromorphicOn.mono_set (Set.subset_univ (closedBall c R))).divisor_ball_support_finite

#print axioms divisor_LFunction_ball_support_finite
