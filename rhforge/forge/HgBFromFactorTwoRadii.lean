import NormProdHalfBall
import NormFactorLe
import NormLeClosedBallOfSphere
import ReLogDerivGeNormalizedAt
import Mathlib.NumberTheory.LSeries.RiemannZeta

/-!
# §1378 — `hgB` from the half-disk factorization, correct two-radii domains

Corrects §1377's domain: `g = ζ/P_half` (only the half-disk zeros removed) is nonvanishing on the *inner* disk
`ball c (R/2)`, not on `ball c R` (it keeps the annulus zeros).  So the ζ-bound and max-modulus use the outer
radius `R`, while the normalized Borel–Carathéodory (§1369) is applied on the inner radius `R/2`, evaluating
`logDeriv g` at the centre `c`.  This is the genuine Titchmarsh §3.9 two-radii construction:

  `-(4 (log M' - log‖g c‖)/(R/2)) ≤ Re (logDeriv g c)`,   `M' = Mζ/(R/2)^{Σ m ρ}`.

Chain: §1376 (`‖P_half‖ ≥ (R/2)^n` on the `R`-sphere) → §1373 (`‖g‖ ≤ M'` on the `R`-sphere) → §1370 (max-modulus
to `closedBall c R ⊇ ball c (R/2)`) → §1369 (normalized BC on radius `R/2`).  = §1355's `hgB` at `c = σ+iγ`.
EXIT0 / std-axioms.
-/

open Complex Metric

theorem re_logDeriv_ge_of_halfBall_factor' {c : ℂ} {R Mζ : ℝ} {g P : ℂ → ℂ} {S : Finset ℂ}
    {m : ℂ → ℕ} (hR : 0 < R) (hPdef : P = fun s => ∏ ρ ∈ S, (s - ρ) ^ (m ρ))
    (hS : ∀ ρ ∈ S, ρ ∈ closedBall c (R / 2))
    (hfac : ∀ s ∈ sphere c R, riemannZeta s = P s * g s)
    (hζsphere : ∀ s ∈ sphere c R, ‖riemannZeta s‖ ≤ Mζ)
    (hgdiff : DiffContOnCl ℂ g (ball c R)) (hgnz : ∀ z ∈ ball c (R / 2), g z ≠ 0)
    (hgc_lt : ‖g c‖ < Mζ / (R / 2) ^ (∑ ρ ∈ S, m ρ)) :
    -(4 * (Real.log (Mζ / (R / 2) ^ (∑ ρ ∈ S, m ρ)) - Real.log ‖g c‖) / (R / 2))
      ≤ (logDeriv g c).re := by
  set M' : ℝ := Mζ / (R / 2) ^ (∑ ρ ∈ S, m ρ) with hM'
  have hR2 : 0 < R / 2 := by linarith
  have hR2pos : 0 < (R / 2) ^ (∑ ρ ∈ S, m ρ) := by positivity
  have hgsphere : ∀ s ∈ sphere c R, ‖g s‖ ≤ M' := by
    intro s hs
    have hPlow : (R / 2) ^ (∑ ρ ∈ S, m ρ) ≤ ‖P s‖ := by
      rw [hPdef]; exact norm_prod_sub_half_ge hR.le m hs hS
    exact norm_factor_le (hfac s hs) hR2pos hPlow (hζsphere s hs)
  have hgdisk : ∀ z ∈ closedBall c R, ‖g z‖ ≤ M' :=
    fun z hz => norm_le_on_closedBall_of_sphere hR hgdiff hgsphere hz
  have hsub : ball c (R / 2) ⊆ ball c R := Metric.ball_subset_ball (by linarith)
  have hgball : ∀ z ∈ ball c (R / 2), ‖g z‖ ≤ M' :=
    fun z hz => hgdisk z (ball_subset_closedBall (hsub hz))
  have hgdon : DifferentiableOn ℂ g (ball c (R / 2)) := hgdiff.differentiableOn.mono hsub
  exact re_logDeriv_ge_of_bounded_normalized_at hR2 hgdon hgnz hgball hgc_lt

#print axioms re_logDeriv_ge_of_halfBall_factor'
