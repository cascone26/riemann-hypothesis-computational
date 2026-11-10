import NormZetaLeCenter
import LogBoundCollapse

/-! [S3020 SHARP COPY: Borel–Carathéodory constant 4 → 2 via SharpBC (S3018); names suffixed _sh]

# §1419 — explicit BC bound `BOUND ≤ 4·log(Mζ/‖ζ c‖)/(R/2)`

§1402's Borel–Carathéodory bound `BOUND = 4(log(Mζ/(R/2)^N) - log‖g_full c‖)/(R/2)` still mentions the opaque
factor `g_full c`.  Feeding §1417 (`‖ζ c‖ ≤ (R/2)^N·‖g_full c‖`) through §1418 collapses the `N`-term, giving the
**explicit** `4·log(Mζ/‖ζ c‖)/(R/2)` — depending only on `Mζ`, `R`, and the (known-function) value `‖ζ c‖`.
This is what lets §1411's `BOUND ≤ B` be stated without any reference to `g₀`.  EXIT0 / std-axioms.
-/

open Metric Classical Filter

theorem bound_le_explicit_sh {c : ℂ} {R' R Mζ : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hball : (1 : ℂ) ∉ ball c R') {g : ℂ → ℂ} (hg : AnalyticOnNhd ℂ g (ball c R'))
    (hcod : riemannZeta =ᶠ[codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u)) • g)
    (hMζ : 0 < Mζ)
    (hGpos : 0 < ‖(∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
        (fun u => ¬ u ∈ closedBall c (R / 2)),
        (c - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g c‖)
    (hZpos : 0 < ‖riemannZeta c‖) :
    2 * (Real.log (Mζ / (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
          (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat))
        - Real.log ‖(∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
            (fun u => ¬ u ∈ closedBall c (R / 2)),
            (c - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g c‖) / (R / 2)
      ≤ 2 * Real.log (Mζ / ‖riemannZeta c‖) / (R / 2) := by
  have hPpos : (0 : ℝ) < (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
      (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) :=
    pow_pos (by linarith) _
  have hlog := log_bound_collapse hMζ hPpos hGpos hZpos (norm_zeta_le_center hR hRR' hball hg hcod)
  have hRhalf : (0 : ℝ) < R / 2 := by linarith
  gcongr

#print axioms bound_le_explicit_sh
