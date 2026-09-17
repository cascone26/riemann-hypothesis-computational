import SmulFactorEqSplitLFunction
import EqOfCodiscreteContinuous
import SplitFunAnalytic
import Mathlib.NumberTheory.LSeries.DirichletContinuation

/-!
# §1505 — centre-value split factorization for `LFunction χ` at ANY ball point (`χ ≠ 1`)

L-analogue of ζ §1416 (cleaner: entire, no `hball`).  §1481 gave the split on the sphere; the same proof works at
any `z ∈ ball c R'`, in particular the centre `c`.  Supplies `LFunction χ c = P_half(c)·g_full(c)` for the
`‖g_full c‖` lower bound.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Classical Filter

theorem LFunction_eq_split_at {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c : ℂ} {R' : ℝ} {g : ℂ → ℂ} (hg : AnalyticOnNhd ℂ g (ball c R'))
    (hcod : LFunction χ =ᶠ[codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u)) • g)
    (p : ℂ → Prop) {z : ℂ} (hz : z ∈ ball c R') :
    LFunction χ z
      = (∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter p,
          (z - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat)
        * ((∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter (fun u => ¬ p u),
            (z - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g z) := by
  rw [smul_factor_eq_split_LFunction hχ p] at hcod
  refine eq_of_codiscreteWithin_of_continuousAt isOpen_ball hz hcod ?_ ?_
  · exact (differentiable_LFunction hχ z).continuousAt
  · exact (split_fun_analyticOnNhd
      (fun u => (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) hg z hz).continuousAt

#print axioms LFunction_eq_split_at
