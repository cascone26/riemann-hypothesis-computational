import HgBConcreteZpow
import ZetaHalfDiskDecompOfFactor
import HgrowthFromDecomp
import DivisorNonnegAll
import DivisorFilterReLe

/-!
# §1404 — concrete `hgrowth` for ζ: `Re(-logDeriv ζ c') ≤ B` at the growth center `c' = σ+2iγ`

Combining the half-disk decomposition §1400 (`logDeriv ζ c' = ∑ divisor/(c'-ρ) + logDeriv g_full c'`) with the
two-radii BC bound §1402 (`-B ≤ Re(logDeriv g_full c')`) through the drop-zeros bound §1399 (using §1387 `divisor ≥
0` and §1403 `Re ρ ≤ Re c'`), the growth term at the second disk center `c'` obeys `Re(-logDeriv ζ c') ≤ B`.  This
is §1336's `hgrowth` with `C = B ≍ log|t|`.  EXIT0 / std-axioms.
-/

open Complex Metric Classical Filter

theorem hgrowth_riemannZeta_concrete {c' : ℂ} {R' R Mζ : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hc're : 1 < c'.re) (hball : (1 : ℂ) ∉ ball c' R') {g₀ : ℂ → ℂ}
    (hga : AnalyticOnNhd ℂ g₀ (ball c' R')) (hgn : ∀ u : ball c' R', g₀ u ≠ 0)
    (hfac_cod : riemannZeta =ᶠ[codiscreteWithin (ball c' R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor riemannZeta (ball c' R') u)) • g₀)
    (hζsphere : ∀ s ∈ sphere c' R, ‖riemannZeta s‖ ≤ Mζ)
    (hgc_lt : ‖(∏ u ∈ (divisor_riemannZeta_ball_support_finite c' R').toFinset.filter
        (fun u => ¬ u ∈ closedBall c' (R / 2)),
        (c' - u) ^ (MeromorphicOn.divisor riemannZeta (ball c' R') u).toNat) * g₀ c'‖
      < Mζ / (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c' R').toFinset.filter
        (· ∈ closedBall c' (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c' R') u).toNat)) :
    (-logDeriv riemannZeta c').re
      ≤ 4 * (Real.log (Mζ / (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c' R').toFinset.filter
          (· ∈ closedBall c' (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c' R') u).toNat))
          - Real.log ‖(∏ u ∈ (divisor_riemannZeta_ball_support_finite c' R').toFinset.filter
            (fun u => ¬ u ∈ closedBall c' (R / 2)),
            (c' - u) ^ (MeromorphicOn.divisor riemannZeta (ball c' R') u).toNat) * g₀ c'‖) / (R / 2) := by
  have hc'ball : c' ∈ ball c' R' := mem_ball_self (by linarith)
  have hdecomp := logDeriv_riemannZeta_halfdisk_decomp_of_factor hc'ball hc're hga hgn hfac_cod
    (fun u => u ∈ closedBall c' (R / 2))
  have hgB := hgB_riemannZeta_concrete_zpow hR hRR' hball hga hgn hfac_cod hζsphere hgc_lt
  exact hgrowth_from_decomp hdecomp (fun ρ _ => divisor_riemannZeta_nonneg_all hball ρ)
    (fun ρ hρ => divisor_filter_support_re_le (le_of_lt hc're) _ hρ) hgB

#print axioms hgrowth_riemannZeta_concrete
