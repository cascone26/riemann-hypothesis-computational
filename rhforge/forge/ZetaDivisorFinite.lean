import Mathlib.NumberTheory.LSeries.RiemannZeta
import Mathlib.Analysis.Meromorphic.Divisor
import ExplicitFormula

/-!
# §1331 — extract_zeros_poles prereq for ζ: `divisor ζ (ball c R)` has finite support

The second ζ-specific hypothesis of `MeromorphicOn.extract_zeros_poles`.  `ζ` is meromorphic on `closedBall c R`
(restriction of `meromorphicOn_riemannZeta`), which is compact (`ℂ` is a `ProperSpace`), so Mathlib's
`divisor_ball_support_finite` gives that ζ has finitely many zeros and poles on the open ball — the divisor has
finite support.  With §1330 (`order ≠ ⊤`) this discharges BOTH ζ-specific hypotheses of the divisor factorization
`ζ = P·g` (`extract_zeros_poles`).  EXIT0 / std-axioms.
-/

open Complex Metric

theorem divisor_riemannZeta_ball_support_finite (c : ℂ) (R : ℝ) :
    (MeromorphicOn.divisor riemannZeta (ball c R)).support.Finite :=
  (meromorphicOn_riemannZeta.mono_set (Set.subset_univ (closedBall c R))).divisor_ball_support_finite

#print axioms divisor_riemannZeta_ball_support_finite
