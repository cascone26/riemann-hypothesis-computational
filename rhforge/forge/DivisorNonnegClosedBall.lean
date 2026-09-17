import AnalyticAtZeta
import ExplicitFormula
import Mathlib.Analysis.Meromorphic.Divisor
import Mathlib.Analysis.Meromorphic.Order

/-!
# §1365 — `0 ≤ divisor ζ (closedBall c R) u` on a pole-free closed ball

The closed-ball companion of §1350 (which was for the open ball).  On `closedBall c |R|` avoiding the pole
`s = 1`, ζ is analytic at every point, so its meromorphic order is `≥ 0` and its divisor is `≥ 0`.  Extracted
as a standalone lemma (previously inlined in §1364) for reuse across the Jensen/zero-count sub-arc.
EXIT0 / std-axioms.
-/

open Complex Metric

theorem divisor_riemannZeta_nonneg_closedBall {c : ℂ} {R : ℝ}
    (hball : (1 : ℂ) ∉ closedBall c |R|) {u : ℂ} (hu : u ∈ closedBall c |R|) :
    0 ≤ MeromorphicOn.divisor riemannZeta (closedBall c |R|) u := by
  have hmero : MeromorphicOn riemannZeta (closedBall c |R|) :=
    meromorphicOn_riemannZeta.mono_set (Set.subset_univ _)
  have hu1 : u ≠ 1 := fun h => hball (h ▸ hu)
  have hord : 0 ≤ meromorphicOrderAt riemannZeta u :=
    (analyticAt_riemannZeta hu1).meromorphicOrderAt_nonneg
  rw [MeromorphicOn.divisor_apply hmero hu]
  rcases eq_or_ne (meromorphicOrderAt riemannZeta u) ⊤ with h | h
  · rw [h]; simp
  · lift meromorphicOrderAt riemannZeta u to ℤ using h with n; simpa using hord

#print axioms divisor_riemannZeta_nonneg_closedBall
