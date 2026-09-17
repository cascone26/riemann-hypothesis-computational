import NormGFullCenter

/-!
# §1414 — discharge §1402's `hgc_lt` from the sphere bound + a strict scalar slack

§1413 gives `‖g_full c‖ ≤ Mζ/(R/2)^N` (non-strict).  §1402 needs the STRICT `hgc_lt`
`‖g_full c‖ < Mζ_h/(R/2)^N`.  Taking the `hgc_lt` constant `Mζ_h` strictly above the actual sphere bound
`Mζ` closes the gap: `‖g_full c‖ ≤ Mζ/(R/2)^N < Mζ_h/(R/2)^N`.  This reduces `hgc_lt` to a sphere bound
(§1398/§1130) plus a scalar slack `Mζ < Mζ_h` — pure fixed-point bookkeeping.  EXIT0 / std-axioms.
-/

open Complex Metric Classical

theorem hgc_lt_of_sphere {c : ℂ} {R' R Mζ Mζ_h : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hball : (1 : ℂ) ∉ ball c R') {g : ℂ → ℂ} (hg : AnalyticOnNhd ℂ g (ball c R'))
    (hcod : riemannZeta =ᶠ[Filter.codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u)) • g)
    (hζsphere : ∀ z ∈ sphere c R, ‖riemannZeta z‖ ≤ Mζ) (hlt : Mζ < Mζ_h) :
    ‖(∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
        (fun u => ¬ u ∈ closedBall c (R / 2)),
        (c - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g c‖
      < Mζ_h / (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
          (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) := by
  have hden : (0 : ℝ) < (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
      (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) :=
    pow_pos (by linarith) _
  calc ‖(∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
          (fun u => ¬ u ∈ closedBall c (R / 2)),
          (c - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g c‖
      ≤ Mζ / (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
          (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) :=
        norm_gfull_le_center hR hRR' hball hg hcod hζsphere
    _ < Mζ_h / (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
          (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) := by
        rw [div_eq_mul_inv, div_eq_mul_inv]
        exact mul_lt_mul_of_pos_right hlt (inv_pos.mpr hden)

#print axioms hgc_lt_of_sphere
