import SmulFactorEqSplit
import EqOfCodiscreteContinuous
import SplitFunAnalytic
import AnalyticAtZeta

/-!
# §1391 — `ζ z = P_half(z)·(P_rest(z)·g z)` on the sphere (§1378's `hfac`)

Rewriting §1343's codiscrete factorization into split form (§1388), then upgrading to pointwise on the sphere via
§1385 — using that both `ζ` (§1345) and the split function (§1390, analytic ⟹ continuous) are continuous at every
`z ∈ sphere c R ⊆ ball c R'` (`R < R'`).  This is exactly the `hfac` hypothesis the two-radii capstone §1378 needs.
EXIT0 / std-axioms.
-/

open Metric Classical Filter

theorem zeta_eq_split_on_sphere {c : ℂ} {R' R : ℝ} (hRR' : R < R') (hball : (1 : ℂ) ∉ ball c R')
    {g : ℂ → ℂ} (hg : AnalyticOnNhd ℂ g (ball c R'))
    (hcod : riemannZeta =ᶠ[codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u)) • g)
    (p : ℂ → Prop) {z : ℂ} (hz : z ∈ sphere c R) :
    riemannZeta z
      = (∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter p,
          (z - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat)
        * ((∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter (fun u => ¬ p u),
            (z - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g z) := by
  rw [smul_factor_eq_split hball p] at hcod
  have hzball : z ∈ ball c R' := by
    rw [mem_sphere] at hz; rw [mem_ball, hz]; exact hRR'
  refine eq_of_codiscreteWithin_of_continuousAt isOpen_ball hzball hcod ?_ ?_
  · have hz1 : z ≠ 1 := fun h => hball (h ▸ hzball)
    exact (analyticAt_riemannZeta hz1).continuousAt
  · exact (split_fun_analyticOnNhd
      (fun u => (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) hg z hzball).continuousAt

#print axioms zeta_eq_split_on_sphere
