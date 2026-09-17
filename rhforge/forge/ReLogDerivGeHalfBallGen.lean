import NormProdHalfBall
import NormFactorLe
import NormLeClosedBallOfSphere
import ReLogDerivGeNormalizedAt

/-!
# §1471 — generic half-disk Borel–Carathéodory factor bound (`riemannZeta` → arbitrary `F`)

The proven ζ brick §1378 `re_logDeriv_ge_of_halfBall_factor'` uses `riemannZeta` only through its two
factorization hypotheses (`F = P·g` on the sphere, `‖F‖ ≤ Mζ` on the sphere) — the analytic engine
`re_logDeriv_ge_of_bounded_normalized_at` (§1369) is entirely generic in `g`.  Generalizing over the outer
function `F : ℂ → ℂ` therefore costs nothing and yields the *shared* two-radii Borel–Carathéodory bound

  `-(4 (log M' - log‖g c‖)/(R/2)) ≤ Re (logDeriv g c)`,   `M' = Mζ / (R/2)^{Σ m ρ}`,

now usable for **both** `F = riemannZeta` (recovering §1378) and `F = LFunction χ` — the exact brick the
Dirichlet-L zero-free region needs to discharge the `Q`/`R` terms of §1470.  EXIT0 / std-axioms.
-/

open Complex Metric

theorem re_logDeriv_ge_of_halfBall_factor_gen {F : ℂ → ℂ} {c : ℂ} {R Mζ : ℝ} {g P : ℂ → ℂ}
    {S : Finset ℂ} {m : ℂ → ℕ} (hR : 0 < R) (hPdef : P = fun s => ∏ ρ ∈ S, (s - ρ) ^ (m ρ))
    (hS : ∀ ρ ∈ S, ρ ∈ closedBall c (R / 2))
    (hfac : ∀ s ∈ sphere c R, F s = P s * g s)
    (hFsphere : ∀ s ∈ sphere c R, ‖F s‖ ≤ Mζ)
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
    exact norm_factor_le (hfac s hs) hR2pos hPlow (hFsphere s hs)
  have hgdisk : ∀ z ∈ closedBall c R, ‖g z‖ ≤ M' :=
    fun z hz => norm_le_on_closedBall_of_sphere hR hgdiff hgsphere hz
  have hsub : ball c (R / 2) ⊆ ball c R := Metric.ball_subset_ball (by linarith)
  have hgball : ∀ z ∈ ball c (R / 2), ‖g z‖ ≤ M' :=
    fun z hz => hgdisk z (ball_subset_closedBall (hsub hz))
  have hgdon : DifferentiableOn ℂ g (ball c (R / 2)) := hgdiff.differentiableOn.mono hsub
  exact re_logDeriv_ge_of_bounded_normalized_at hR2 hgdon hgnz hgball hgc_lt

#print axioms re_logDeriv_ge_of_halfBall_factor_gen
