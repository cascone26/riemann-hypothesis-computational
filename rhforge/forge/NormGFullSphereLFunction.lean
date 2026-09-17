import LFunctionEqSplitSphere
import NormProdHalfBall
import NormFactorLe

/-!
# §1502 — sphere bound on `g_full` for `LFunction χ`: `‖g_full z‖ ≤ Mζ/(R/2)^N` (`z ∈ sphere c R`)

L-analogue of ζ §1412 (cleaner: entire, no `hball`).  On the outer sphere, `LFunction χ = P_half · g_full` (§1481),
the half-disk product has `‖P_half z‖ ≥ (R/2)^N` (generic §1376), and `‖LFunction χ z‖ ≤ Mζ`, so generic §1373
gives `‖g_full z‖ ≤ Mζ/(R/2)^N`.  Sphere-level input for the max-modulus step that discharges the centre `hgc_lt`.
EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Classical

theorem norm_gfull_le_sphere_LFunction {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c z : ℂ} {R' R Mζ : ℝ} (hR : 0 < R) (hRR' : R < R')
    {g : ℂ → ℂ} (hg : AnalyticOnNhd ℂ g (ball c R'))
    (hcod : LFunction χ =ᶠ[Filter.codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u)) • g)
    (hζ : ‖LFunction χ z‖ ≤ Mζ) (hz : z ∈ sphere c R) :
    ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
        (fun u => ¬ u ∈ closedBall c (R / 2)),
        (z - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g z‖
      ≤ Mζ / (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
          (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) := by
  have hsplit := LFunction_eq_split_on_sphere hχ hRR' hg hcod (fun u => u ∈ closedBall c (R / 2)) hz
  have hPlow := norm_prod_sub_half_ge
    (S := (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
      (· ∈ closedBall c (R / 2))) hR.le
    (fun u => (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) hz
    (fun ρ hρ => (Finset.mem_filter.mp hρ).2)
  have hPpos : (0 : ℝ) < (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
      (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) :=
    pow_pos (by linarith) _
  exact norm_factor_le hsplit hPpos hPlow hζ

#print axioms norm_gfull_le_sphere_LFunction
