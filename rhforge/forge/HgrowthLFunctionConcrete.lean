import LFunctionHalfDiskDecomp
import HgBLFunctionConcreteZpow
import HgrowthFromDecomp
import DivisorLFunctionNonneg
import DivisorFilterReLeLFunction

/-!
# §1491 — concrete `hgrowth` for `LFunction χ`: `Re(-logDeriv (LFunction χ) c') ≤ B` at `c' = σ+2iγ`

L-analogue of ζ §1404 (cleaner: entire, no `hball`).  Combining the half-disk decomposition §1486
(`logDeriv (LFunction χ) c' = ∑ divisor/(c'-ρ) + logDeriv g_full c'`) with the ℤ-exponent two-radii BC bound §1488
(`-B ≤ Re(logDeriv g_full c')`) through the generic drop-zeros bound §1399 (using §1479 `divisor ≥ 0` and §1490
`Re ρ ≤ Re c'`), the growth term at the second-disk center `c'` obeys `Re(-logDeriv (LFunction χ) c') ≤ B`, with
`B ≍ log|t|`.  For `χ² ≠ 1` this discharges the `R` term of §1470 directly; the same brick applied at `c' = σ+iγ`
gives the drop-all bound feeding the `Q` (keep-one-zero) step.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Classical Filter

theorem hgrowth_LFunction_concrete {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c' : ℂ} {R' R Mζ : ℝ} (hR : 0 < R) (hRR' : R < R') (hc're : 1 < c'.re)
    {g₀ : ℂ → ℂ} (hga : AnalyticOnNhd ℂ g₀ (ball c' R')) (hgn : ∀ u : ball c' R', g₀ u ≠ 0)
    (hfac_cod : LFunction χ =ᶠ[codiscreteWithin (ball c' R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c' R') u)) • g₀)
    (hLsphere : ∀ s ∈ sphere c' R, ‖LFunction χ s‖ ≤ Mζ)
    (hgc_lt : ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ c' R').toFinset.filter
        (fun u => ¬ u ∈ closedBall c' (R / 2)),
        (c' - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c' R') u).toNat) * g₀ c'‖
      < Mζ / (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c' R').toFinset.filter
        (· ∈ closedBall c' (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c' R') u).toNat)) :
    (-logDeriv (LFunction χ) c').re
      ≤ 4 * (Real.log (Mζ / (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c' R').toFinset.filter
          (· ∈ closedBall c' (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c' R') u).toNat))
          - Real.log ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ c' R').toFinset.filter
            (fun u => ¬ u ∈ closedBall c' (R / 2)),
            (c' - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c' R') u).toNat) * g₀ c'‖) / (R / 2) := by
  have hc'ball : c' ∈ ball c' R' := mem_ball_self (by linarith)
  have hdecomp := logDeriv_LFunction_halfdisk_decomp_of_factor hχ hc'ball hc're hga hgn hfac_cod
    (fun u => u ∈ closedBall c' (R / 2))
  have hgB := hgB_LFunction_concrete_zpow hχ hR hRR' hga hgn hfac_cod hLsphere hgc_lt
  exact hgrowth_from_decomp hdecomp (fun ρ _ => divisor_LFunction_nonneg_all hχ ρ)
    (fun ρ hρ => divisor_filter_support_re_le_LFunction hχ (le_of_lt hc're) _ hρ) hgB

#print axioms hgrowth_LFunction_concrete
