import LFunctionEqSplitAt
import NormProdSubCenterLe

/-!
# §1506 — centre value bound `‖LFunction χ c‖ ≤ (R/2)^N · ‖g_full c‖` (`χ ≠ 1`)

L-analogue of ζ §1417.  Evaluate the centre split (§1505) `LFunction χ c = P_half(c)·g_full(c)`, take norms, bound
`‖P_half c‖ ≤ (R/2)^N` (generic §1415).  Gives the `‖g_full c‖` lower bound that collapses the BC `BOUND`.
EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Classical Filter

theorem norm_LFunction_le_center {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c : ℂ} {R' R : ℝ} (hR : 0 < R) (hRR' : R < R')
    {g : ℂ → ℂ} (hg : AnalyticOnNhd ℂ g (ball c R'))
    (hcod : LFunction χ =ᶠ[codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u)) • g) :
    ‖LFunction χ c‖
      ≤ (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
          (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat)
        * ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
            (fun u => ¬ u ∈ closedBall c (R / 2)),
            (c - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g c‖ := by
  have hc : c ∈ ball c R' := mem_ball_self (lt_trans hR hRR')
  rw [LFunction_eq_split_at hχ hg hcod (· ∈ closedBall c (R / 2)) hc, norm_mul]
  exact mul_le_mul_of_nonneg_right
    (norm_prod_sub_center_le (by linarith) _ (fun ρ hρ => (Finset.mem_filter.mp hρ).2))
    (norm_nonneg _)

#print axioms norm_LFunction_le_center
