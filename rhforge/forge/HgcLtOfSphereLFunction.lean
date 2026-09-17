import NormGFullCenterLFunction

/-!
# §1504 — discharge the centre `hgc_lt` for `LFunction χ` from the sphere bound + strict slack

L-analogue of ζ §1414 (cleaner: entire, no `hball`).  §1503 gives `‖g_full c‖ ≤ Mζ/(R/2)^N` (non-strict); the
strict `hgc_lt` needed by §1493/§1488 follows by taking the `hgc_lt` constant `Mζ_h` strictly above the sphere
bound `Mζ`: `‖g_full c‖ ≤ Mζ/(R/2)^N < Mζ_h/(R/2)^N`.  Reduces `hgc_lt` to a sphere bound + a scalar slack
`Mζ < Mζ_h` (pure fixed-point bookkeeping).  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Classical

theorem hgc_lt_of_sphere_LFunction {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c : ℂ} {R' R Mζ Mζ_h : ℝ} (hR : 0 < R) (hRR' : R < R')
    {g : ℂ → ℂ} (hg : AnalyticOnNhd ℂ g (ball c R'))
    (hcod : LFunction χ =ᶠ[Filter.codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u)) • g)
    (hζsphere : ∀ z ∈ sphere c R, ‖LFunction χ z‖ ≤ Mζ) (hlt : Mζ < Mζ_h) :
    ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
        (fun u => ¬ u ∈ closedBall c (R / 2)),
        (c - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g c‖
      < Mζ_h / (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
          (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) := by
  have hden : (0 : ℝ) < (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
      (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) :=
    pow_pos (by linarith) _
  calc ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
          (fun u => ¬ u ∈ closedBall c (R / 2)),
          (c - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g c‖
      ≤ Mζ / (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
          (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) :=
        norm_gfull_le_center_LFunction hχ hR hRR' hg hcod hζsphere
    _ < Mζ_h / (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
          (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) := by
        rw [div_eq_mul_inv, div_eq_mul_inv]
        exact mul_lt_mul_of_pos_right hlt (inv_pos.mpr hden)

#print axioms hgc_lt_of_sphere_LFunction
