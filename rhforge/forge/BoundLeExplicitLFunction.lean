import NormLFunctionLeCenter
import LogBoundCollapse

/-!
# §1507 — explicit BC bound for `LFunction χ`: `BOUND ≤ 4·log(Mζ/‖L c‖)/(R/2)`

L-analogue of ζ §1419.  The BC `BOUND = 4(log(Mζ/(R/2)^N) - log‖g_full c‖)/(R/2)` still mentions the opaque
`g_full c`.  Feeding §1506 (`‖L c‖ ≤ (R/2)^N·‖g_full c‖`) through the generic collapse §1418 gives the explicit
`4·log(Mζ/‖L c‖)/(R/2)` — depending only on `Mζ`, `R`, and the value `‖LFunction χ c‖`.  This lets `BOUND ≤ B` be
stated with no reference to `g₀`.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Classical Filter

theorem bound_le_explicit_LFunction {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c : ℂ} {R' R Mζ : ℝ} (hR : 0 < R) (hRR' : R < R')
    {g : ℂ → ℂ} (hg : AnalyticOnNhd ℂ g (ball c R'))
    (hcod : LFunction χ =ᶠ[codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u)) • g)
    (hMζ : 0 < Mζ)
    (hGpos : 0 < ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
        (fun u => ¬ u ∈ closedBall c (R / 2)),
        (c - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g c‖)
    (hZpos : 0 < ‖LFunction χ c‖) :
    4 * (Real.log (Mζ / (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
          (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat))
        - Real.log ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
            (fun u => ¬ u ∈ closedBall c (R / 2)),
            (c - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g c‖) / (R / 2)
      ≤ 4 * Real.log (Mζ / ‖LFunction χ c‖) / (R / 2) := by
  have hPpos : (0 : ℝ) < (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
      (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) :=
    pow_pos (by linarith) _
  have hlog := log_bound_collapse hMζ hPpos hGpos hZpos
    (norm_LFunction_le_center hχ hR hRR' hg hcod)
  have hRhalf : (0 : ℝ) < R / 2 := by linarith
  gcongr

#print axioms bound_le_explicit_LFunction
