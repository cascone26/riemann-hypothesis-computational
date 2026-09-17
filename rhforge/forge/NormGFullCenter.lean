import NormLeClosedBallOfSphere
import NormGFullSphere
import GFullAnalytic
import AnalyticDiffContOnCl

/-!
# §1413 — max-modulus centre bound: `‖g_full c‖ ≤ Mζ/(R/2)^N`

`g_full = (∏_{¬half}(·-u)^toNat)·g₀` is analytic on `ball c R'` (§1383), hence `DiffContOnCl` on the smaller
`ball c R` (§1384, `R<R'`); its sphere bound is `Mζ/(R/2)^N` (§1412).  The maximum-modulus principle (§1370)
propagates that bound to the whole closed disk, in particular to the centre `c`.  This is the value bound whose
strict slack discharges §1402's `hgc_lt`.  EXIT0 / std-axioms.
-/

open Complex Metric Classical

theorem norm_gfull_le_center {c : ℂ} {R' R Mζ : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hball : (1 : ℂ) ∉ ball c R') {g : ℂ → ℂ} (hg : AnalyticOnNhd ℂ g (ball c R'))
    (hcod : riemannZeta =ᶠ[Filter.codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u)) • g)
    (hζsphere : ∀ z ∈ sphere c R, ‖riemannZeta z‖ ≤ Mζ) :
    ‖(∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
        (fun u => ¬ u ∈ closedBall c (R / 2)),
        (c - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g c‖
      ≤ Mζ / (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
          (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) := by
  have han := prod_rest_mul_analyticOnNhd
    (T := (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
      (fun u => ¬ u ∈ closedBall c (R / 2)))
    (fun u => (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) hg
  have hdc := analyticOnNhd_diffContOnCl hR hRR' han
  exact norm_le_on_closedBall_of_sphere hR hdc
    (fun z hz => norm_gfull_le_sphere hR hRR' hball hg hcod (hζsphere z hz) hz)
    (mem_closedBall_self hR.le)

#print axioms norm_gfull_le_center
