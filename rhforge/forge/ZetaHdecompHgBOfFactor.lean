import ZetaHalfDiskDecompOfFactor
import HgBConcreteZpow

/-!
# §1410 — coordinated `hdecomp ∧ hgB` from ONE §1343 factorization (§1409's two analytic slots share `g`)

From a single §1343 factorization of ζ on `ball c R'` (plus the sphere bound `Mζ` on `sphere c R`, `R<R'`, and the
strict factor bound `hgc_lt`), produces the **same** factor value `g = logDeriv g_full c` serving BOTH of §1409's
analytic slots at once:
* `hdecomp` — the half-disk logarithmic-derivative decomposition (§1400 at radius `R'`, `p = (·∈closedBall c (R/2))`);
* `hgB` — the Borel–Carathéodory lower bound `-bound ≤ g.re` (§1402).

This is the coordination content — the two bounds must be about the *identical* `g_full` built from the *identical*
`g₀`.  `c` is kept abstract (only `1 < c.re` needed).  EXIT0 / std-axioms.
-/

open Complex Metric Filter Classical

theorem zeta_hdecomp_and_hgB_of_factor {c : ℂ} {R' R Mζ : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hcre : 1 < c.re) (hball : (1 : ℂ) ∉ ball c R') {g₀ : ℂ → ℂ}
    (hga : AnalyticOnNhd ℂ g₀ (ball c R')) (hgn : ∀ u : ball c R', g₀ u ≠ 0)
    (hfac_cod : riemannZeta =ᶠ[codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u)) • g₀)
    (hζsphere : ∀ s ∈ sphere c R, ‖riemannZeta s‖ ≤ Mζ)
    (hgc_lt : ‖(∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
        (fun u => ¬ u ∈ closedBall c (R / 2)),
        (c - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g₀ c‖
      < Mζ / (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
        (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat)) :
    ∃ g : ℂ,
      logDeriv riemannZeta c
        = (∑ ρ ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
            (fun u => u ∈ closedBall c (R / 2)),
            (MeromorphicOn.divisor riemannZeta (ball c R') ρ : ℂ) / (c - ρ)) + g
      ∧ -(4 * (Real.log (Mζ / (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
          (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat))
          - Real.log ‖(∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
            (fun u => ¬ u ∈ closedBall c (R / 2)),
            (c - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g₀ c‖) / (R / 2)) ≤ g.re := by
  have hR' : 0 < R' := lt_trans hR hRR'
  refine ⟨logDeriv (fun s => (∏ ρ ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
      (fun u => ¬ u ∈ closedBall c (R / 2)),
      (s - ρ) ^ (MeromorphicOn.divisor riemannZeta (ball c R') ρ)) * g₀ s) c, ?_, ?_⟩
  · exact logDeriv_riemannZeta_halfdisk_decomp_of_factor (mem_ball_self hR') hcre hga hgn hfac_cod
      (fun u => u ∈ closedBall c (R / 2))
  · exact hgB_riemannZeta_concrete_zpow hR hRR' hball hga hgn hfac_cod hζsphere hgc_lt

#print axioms zeta_hdecomp_and_hgB_of_factor
