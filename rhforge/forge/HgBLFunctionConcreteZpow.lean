import HgBLFunctionConcrete
import GFullZpowEqNatLFunction

/-!
# §1488 — concrete `hgB` about the **ℤ-exponent** `g_full` for `LFunction χ` (matches §1486's decomposition)

L-analogue of ζ §1402.  §1482 bounds `Re(logDeriv g_full c)` for the ℕ-exponent factor.  §1487 says the ℤ- and
ℕ-exponent `g_full` are the same function, so their `logDeriv`s agree — rewriting gives the identical bound about
the **ℤ-exponent** factor, exactly the `g_full` appearing in §1486's half-disk decomposition.  This is the `hgB`
(`-B ≤ Re g`, `g = logDeriv g_full c`) that composes with §1486 to discharge the DLVP `Q`/`R` terms.
EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Classical

theorem hgB_LFunction_concrete_zpow {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c : ℂ} {R' R Mζ : ℝ} (hR : 0 < R) (hRR' : R < R')
    {g₀ : ℂ → ℂ} (hga : AnalyticOnNhd ℂ g₀ (ball c R'))
    (hgn : ∀ u : ball c R', g₀ u ≠ 0)
    (hfac_cod : LFunction χ =ᶠ[Filter.codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u)) • g₀)
    (hLsphere : ∀ s ∈ sphere c R, ‖LFunction χ s‖ ≤ Mζ)
    (hgc_lt : ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
        (fun u => ¬ u ∈ closedBall c (R / 2)),
        (c - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g₀ c‖
      < Mζ / (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
        (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat)) :
    -(4 * (Real.log (Mζ / (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
        (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat))
        - Real.log ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
          (fun u => ¬ u ∈ closedBall c (R / 2)),
          (c - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g₀ c‖) / (R / 2))
      ≤ (logDeriv (fun s => (∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
          (fun u => ¬ u ∈ closedBall c (R / 2)),
          (s - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u)) * g₀ s) c).re := by
  rw [gfull_zpow_eq_nat_LFunction (hχ := hχ)
    (T := (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
      (fun u => ¬ u ∈ closedBall c (R / 2))) (g := g₀)]
  exact hgB_LFunction_concrete hχ hR hRR' hga hgn hfac_cod hLsphere hgc_lt

#print axioms hgB_LFunction_concrete_zpow
