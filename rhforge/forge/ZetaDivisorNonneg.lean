import AnalyticAtZeta
import ExplicitFormula
import Mathlib.Analysis.Meromorphic.Divisor
import Mathlib.Analysis.Meromorphic.Order

/-!
# §1350 — `0 ≤ divisor ζ (ball c R) ρ` on a pole-free disk

On a disk `B(c,R)` avoiding the pole `s = 1`, ζ is analytic at each point (§1345), so its meromorphic order is
nonnegative (`AnalyticAt.meromorphicOrderAt_nonneg`), and `divisor = order.untop₀ ≥ 0`.  Every divisor value is a
nonnegative multiplicity — §1335's `hm : ∀ ρ, 0 ≤ m ρ`.  EXIT0 / std-axioms.
-/

open Complex Metric

theorem divisor_riemannZeta_nonneg {c : ℂ} {R : ℝ} (hball : (1 : ℂ) ∉ ball c R)
    {ρ : ℂ} (hρ : ρ ∈ ball c R) :
    0 ≤ MeromorphicOn.divisor riemannZeta (ball c R) ρ := by
  have hmero : MeromorphicOn riemannZeta (ball c R) :=
    meromorphicOn_riemannZeta.mono_set (Set.subset_univ _)
  have hρ1 : ρ ≠ 1 := fun h => hball (h ▸ hρ)
  have hord : 0 ≤ meromorphicOrderAt riemannZeta ρ :=
    (analyticAt_riemannZeta hρ1).meromorphicOrderAt_nonneg
  rw [MeromorphicOn.divisor_apply hmero hρ]
  rcases eq_or_ne (meromorphicOrderAt riemannZeta ρ) ⊤ with h | h
  · rw [h]; simp
  · lift meromorphicOrderAt riemannZeta ρ to ℤ using h with n
    simpa using hord

#print axioms divisor_riemannZeta_nonneg
