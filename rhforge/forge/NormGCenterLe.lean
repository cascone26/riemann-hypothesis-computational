import NormLogDerivQuarterBall

set_option maxHeartbeats 1000000

open Complex Metric

/-!
# §1682 — `‖g c‖ ≤ M₀/(3R/4)^Σm`: discharging the master's strict `hgc_lt` with `M`-slack

The center bound for the zero-free factor, extracted standalone from §1679's first half (sphere
division + max modulus). Feeding the master §1680 with `M := M₀ + 1` (any strict slack) turns this
`≤` into the required strict `<`, since `M₀/(3R/4)^Σ < (M₀+1)/(3R/4)^Σ`.
-/

theorem norm_g_center_le_of_quarterBall_factor {c : ℂ} {R M₀ : ℝ} {g P : ℂ → ℂ} {S : Finset ℂ}
    {m : ℂ → ℕ} (hR : 0 < R) (hPdef : P = fun s => ∏ ρ ∈ S, (s - ρ) ^ (m ρ))
    (hS : ∀ ρ ∈ S, ρ ∈ closedBall c (R / 4))
    (hfac : ∀ s ∈ sphere c R, riemannZeta s = P s * g s)
    (hζsphere : ∀ s ∈ sphere c R, ‖riemannZeta s‖ ≤ M₀)
    (hgdiff : DiffContOnCl ℂ g (ball c R)) :
    ‖g c‖ ≤ M₀ / (3 * R / 4) ^ (∑ ρ ∈ S, m ρ) := by
  have h34 : (0 : ℝ) < 3 * R / 4 := by linarith
  have hgsphere : ∀ s ∈ sphere c R, ‖g s‖ ≤ M₀ / (3 * R / 4) ^ (∑ ρ ∈ S, m ρ) := by
    intro s hs
    have hPlow : (3 * R / 4) ^ (∑ ρ ∈ S, m ρ) ≤ ‖P s‖ := by
      rw [hPdef, norm_prod]
      simp_rw [norm_pow]
      rw [← Finset.prod_pow_eq_pow_sum]
      apply Finset.prod_le_prod
      · intro ρ _
        positivity
      · intro ρ hρ
        have hbase : 3 * R / 4 ≤ ‖s - ρ‖ := by
          have hρc : dist ρ c ≤ R / 4 := mem_closedBall.mp (hS ρ hρ)
          have hsc : dist s c = R := mem_sphere.mp hs
          have htri := dist_triangle s ρ c
          rw [← dist_eq_norm]
          linarith
        gcongr
    exact norm_factor_le (hfac s hs) (by positivity) hPlow (hζsphere s hs)
  have hfr : ∀ u ∈ frontier (ball c R), ‖g u‖ ≤ M₀ / (3 * R / 4) ^ (∑ ρ ∈ S, m ρ) := by
    intro u hu
    rw [frontier_ball c hR.ne'] at hu
    exact hgsphere u hu
  exact Complex.norm_le_of_forall_mem_frontier_norm_le isBounded_ball hgdiff hfr
    (subset_closure (mem_ball_self hR))

#print axioms norm_g_center_le_of_quarterBall_factor
