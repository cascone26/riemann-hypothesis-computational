import NormProdHalfBall
import NormFactorLe
import NormLeClosedBallOfSphere
import ReLogDerivGeNormalizedAt
import Mathlib.NumberTheory.LSeries.RiemannZeta

/-!
# §1377 — `hgB` from the half-disk factorization (capstone of the factor-bound chain)

Threads the whole two-radii factor bound end-to-end: given the half-disk factorization `ζ = P_half·g` on the
`R`-sphere (with `P_half = ∏_{ρ∈S}(s-ρ)^{m ρ}`, `S ⊆ closedBall c (R/2)`), a ζ-bound `Mζ` on the sphere, and `g`
holomorphic/nonvanishing/`DiffContOnCl` on `ball c R`, the normalized Borel–Carathéodory bound at the centre gives

  `-(4 (log M' - log‖g c‖)/R) ≤ Re (logDeriv g c)`,   `M' = Mζ / (R/2)^{Σ m ρ}`.

Chain: §1376 (`‖P_half‖ ≥ (R/2)^n` on sphere) → §1373 (`‖g‖ ≤ Mζ/(R/2)^n` on sphere) → §1370 (max-modulus to the
disk) → §1369 (normalized BC).  This is exactly §1355's `hgB` at `c = σ+iγ`, `B = 4(log M' - log‖g c‖)/R ≍ log|t|`.
EXIT0 / std-axioms.
-/

open Complex Metric

theorem re_logDeriv_ge_of_halfBall_factor {c : ℂ} {R Mζ : ℝ} {g P : ℂ → ℂ} {S : Finset ℂ}
    {m : ℂ → ℕ} (hR : 0 < R) (hPdef : P = fun s => ∏ ρ ∈ S, (s - ρ) ^ (m ρ))
    (hS : ∀ ρ ∈ S, ρ ∈ closedBall c (R / 2))
    (hfac : ∀ s ∈ sphere c R, riemannZeta s = P s * g s)
    (hζsphere : ∀ s ∈ sphere c R, ‖riemannZeta s‖ ≤ Mζ)
    (hgdiff : DiffContOnCl ℂ g (ball c R)) (hgnz : ∀ z ∈ ball c R, g z ≠ 0)
    (hgc_lt : ‖g c‖ < Mζ / (R / 2) ^ (∑ ρ ∈ S, m ρ)) :
    -(4 * (Real.log (Mζ / (R / 2) ^ (∑ ρ ∈ S, m ρ)) - Real.log ‖g c‖) / R)
      ≤ (logDeriv g c).re := by
  set M' : ℝ := Mζ / (R / 2) ^ (∑ ρ ∈ S, m ρ) with hM'
  have hR2pos : 0 < (R / 2) ^ (∑ ρ ∈ S, m ρ) := by positivity
  have hgsphere : ∀ s ∈ sphere c R, ‖g s‖ ≤ M' := by
    intro s hs
    have hPlow : (R / 2) ^ (∑ ρ ∈ S, m ρ) ≤ ‖P s‖ := by
      rw [hPdef]; exact norm_prod_sub_half_ge hR.le m hs hS
    exact norm_factor_le (hfac s hs) hR2pos hPlow (hζsphere s hs)
  have hgdisk : ∀ z ∈ closedBall c R, ‖g z‖ ≤ M' :=
    fun z hz => norm_le_on_closedBall_of_sphere hR hgdiff hgsphere hz
  have hgball : ∀ z ∈ ball c R, ‖g z‖ ≤ M' :=
    fun z hz => hgdisk z (ball_subset_closedBall hz)
  exact re_logDeriv_ge_of_bounded_normalized_at hR hgdiff.differentiableOn hgnz hgball hgc_lt

#print axioms re_logDeriv_ge_of_halfBall_factor
