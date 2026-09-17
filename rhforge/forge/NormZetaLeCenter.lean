import ZetaEqSplitAt
import NormProdSubCenterLe

/-!
# §1417 — centre value bound `‖ζ c‖ ≤ (R/2)^N · ‖g_full c‖`

Evaluate the centre factorization (§1416) `ζ c = P_half(c)·g_full(c)`, take norms, and bound the half-disk
product `‖P_half c‖ ≤ (R/2)^N` (§1415).  Multiplicative form (no division), giving the `‖g_full c‖` LOWER bound
that makes §1402's `BOUND` collapse to `4·log(Mζ/‖ζ c‖)/(R/2)`.  EXIT0 / std-axioms.
-/

open Metric Classical Filter

theorem norm_zeta_le_center {c : ℂ} {R' R : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hball : (1 : ℂ) ∉ ball c R') {g : ℂ → ℂ} (hg : AnalyticOnNhd ℂ g (ball c R'))
    (hcod : riemannZeta =ᶠ[codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u)) • g) :
    ‖riemannZeta c‖
      ≤ (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
          (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat)
        * ‖(∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
            (fun u => ¬ u ∈ closedBall c (R / 2)),
            (c - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g c‖ := by
  have hc : c ∈ ball c R' := mem_ball_self (lt_trans hR hRR')
  rw [zeta_eq_split_at hball hg hcod (· ∈ closedBall c (R / 2)) hc, norm_mul]
  exact mul_le_mul_of_nonneg_right
    (norm_prod_sub_center_le (by linarith) _ (fun ρ hρ => (Finset.mem_filter.mp hρ).2))
    (norm_nonneg _)

#print axioms norm_zeta_le_center
