import HgBConcrete
import GFullZpowEqNat

/-!
# §1402 — concrete `hgB` about the **ℤ-exponent** `g_full` (matches §1400's decomposition)

§1401 bounds `Re(logDeriv g_full c)` for the ℕ-exponent factor.  §1394 says the ℤ- and ℕ-exponent `g_full` are the
same function, so their `logDeriv`s agree — rewriting gives the identical bound about the **ℤ-exponent** factor,
which is exactly the `g_full` appearing in §1400's half-disk decomposition.  This is §1336's `hgB` (`-B ≤ Re g`,
`g = logDeriv g_full c`) for the real ζ.  EXIT0 / std-axioms.
-/

open Complex Metric Classical

theorem hgB_riemannZeta_concrete_zpow {c : ℂ} {R' R Mζ : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hball : (1 : ℂ) ∉ ball c R') {g₀ : ℂ → ℂ} (hga : AnalyticOnNhd ℂ g₀ (ball c R'))
    (hgn : ∀ u : ball c R', g₀ u ≠ 0)
    (hfac_cod : riemannZeta =ᶠ[Filter.codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u)) • g₀)
    (hζsphere : ∀ s ∈ sphere c R, ‖riemannZeta s‖ ≤ Mζ)
    (hgc_lt : ‖(∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
        (fun u => ¬ u ∈ closedBall c (R / 2)),
        (c - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g₀ c‖
      < Mζ / (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
        (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat)) :
    -(4 * (Real.log (Mζ / (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
        (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat))
        - Real.log ‖(∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
          (fun u => ¬ u ∈ closedBall c (R / 2)),
          (c - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g₀ c‖) / (R / 2))
      ≤ (logDeriv (fun s => (∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
          (fun u => ¬ u ∈ closedBall c (R / 2)),
          (s - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u)) * g₀ s) c).re := by
  rw [gfull_zpow_eq_nat (T := (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
    (fun u => ¬ u ∈ closedBall c (R / 2))) hball (g := g₀)]
  exact hgB_riemannZeta_concrete hR hRR' hball hga hgn hfac_cod hζsphere hgc_lt

#print axioms hgB_riemannZeta_concrete_zpow
