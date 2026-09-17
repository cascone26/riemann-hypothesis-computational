import ZetaEqSplitSphere
import NormProdHalfBall
import NormFactorLe

/-!
# §1412 — sphere bound on `g_full`: `‖g_full z‖ ≤ Mζ/(R/2)^N` for `z ∈ sphere c R`

On the outer sphere, `ζ = P_half · g_full` (§1391), the half-disk product has `‖P_half z‖ ≥ (R/2)^N`
(§1376, `N = ∑ divisor.toNat` over the half-disk zeros), and `‖ζ z‖ ≤ Mζ`, so §1373 gives
`‖g_full z‖ ≤ Mζ/(R/2)^N`.  This is the sphere-level input for the max-modulus step that will discharge
§1402's `hgc_lt` at the centre.  EXIT0 / std-axioms.
-/

open Complex Metric Classical

theorem norm_gfull_le_sphere {c z : ℂ} {R' R Mζ : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hball : (1 : ℂ) ∉ ball c R') {g : ℂ → ℂ} (hg : AnalyticOnNhd ℂ g (ball c R'))
    (hcod : riemannZeta =ᶠ[Filter.codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u)) • g)
    (hζ : ‖riemannZeta z‖ ≤ Mζ) (hz : z ∈ sphere c R) :
    ‖(∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
        (fun u => ¬ u ∈ closedBall c (R / 2)),
        (z - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g z‖
      ≤ Mζ / (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
          (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) := by
  have hsplit := zeta_eq_split_on_sphere hRR' hball hg hcod (fun u => u ∈ closedBall c (R / 2)) hz
  have hPlow := norm_prod_sub_half_ge (S := (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
      (· ∈ closedBall c (R / 2))) hR.le
    (fun u => (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) hz
    (fun ρ hρ => (Finset.mem_filter.mp hρ).2)
  have hPpos : (0 : ℝ) < (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
      (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) :=
    pow_pos (by linarith) _
  exact norm_factor_le hsplit hPpos hPlow hζ

#print axioms norm_gfull_le_sphere
