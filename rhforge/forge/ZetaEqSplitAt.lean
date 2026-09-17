import SmulFactorEqSplit
import EqOfCodiscreteContinuous
import SplitFunAnalytic
import AnalyticAtZeta

/-!
# §1416 — centre-value split factorization `ζ z = P_p(z)·(P_{¬p}(z)·g z)` at ANY ball point

§1391 gave the split factorization on the sphere; its proof only used `z ∈ ball c R'`, so it generalizes verbatim
to any point of the ball — in particular the centre `c`.  This supplies the centre-value equation
`ζ c = P_half(c)·g_full(c)` needed (with §1415) for the `‖g_full c‖ ≥ ‖ζ c‖/r^N` lower bound.  EXIT0 / std-axioms.
-/

open Metric Classical Filter

theorem zeta_eq_split_at {c : ℂ} {R' : ℝ} (hball : (1 : ℂ) ∉ ball c R') {g : ℂ → ℂ}
    (hg : AnalyticOnNhd ℂ g (ball c R'))
    (hcod : riemannZeta =ᶠ[codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u)) • g)
    (p : ℂ → Prop) {z : ℂ} (hz : z ∈ ball c R') :
    riemannZeta z
      = (∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter p,
          (z - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat)
        * ((∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter (fun u => ¬ p u),
            (z - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g z) := by
  rw [smul_factor_eq_split hball p] at hcod
  refine eq_of_codiscreteWithin_of_continuousAt isOpen_ball hz hcod ?_ ?_
  · exact (analyticAt_riemannZeta (fun h => hball (h ▸ hz))).continuousAt
  · exact (split_fun_analyticOnNhd
      (fun u => (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) hg z hz).continuousAt

#print axioms zeta_eq_split_at
