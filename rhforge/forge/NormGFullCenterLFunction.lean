import NormLeClosedBallOfSphere
import NormGFullSphereLFunction
import GFullAnalytic
import AnalyticDiffContOnCl

/-!
# §1503 — max-modulus centre bound for `LFunction χ`: `‖g_full c‖ ≤ Mζ/(R/2)^N`

L-analogue of ζ §1413 (cleaner: entire, no `hball`).  `g_full = (∏_{¬half}(·-u)^toNat)·g₀` is analytic on
`ball c R'` (generic §1383), hence `DiffContOnCl` on `ball c R` (generic §1384); its sphere bound is `Mζ/(R/2)^N`
(§1502).  The maximum-modulus principle (generic §1370) propagates that bound to the centre `c`.  This is the value
bound whose strict slack discharges the centre `hgc_lt`.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Classical

theorem norm_gfull_le_center_LFunction {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c : ℂ} {R' R Mζ : ℝ} (hR : 0 < R) (hRR' : R < R')
    {g : ℂ → ℂ} (hg : AnalyticOnNhd ℂ g (ball c R'))
    (hcod : LFunction χ =ᶠ[Filter.codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u)) • g)
    (hζsphere : ∀ z ∈ sphere c R, ‖LFunction χ z‖ ≤ Mζ) :
    ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
        (fun u => ¬ u ∈ closedBall c (R / 2)),
        (c - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g c‖
      ≤ Mζ / (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
          (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) := by
  have han := prod_rest_mul_analyticOnNhd
    (T := (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
      (fun u => ¬ u ∈ closedBall c (R / 2)))
    (fun u => (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) hg
  have hdc := analyticOnNhd_diffContOnCl hR hRR' han
  exact norm_le_on_closedBall_of_sphere hR hdc
    (fun z hz => norm_gfull_le_sphere_LFunction hχ hR hRR' hg hcod (hζsphere z hz) hz)
    (mem_closedBall_self hR.le)

#print axioms norm_gfull_le_center_LFunction
