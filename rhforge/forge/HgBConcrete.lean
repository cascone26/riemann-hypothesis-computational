import HgBFromFactorTwoRadii
import ZetaEqSplitSphere
import GFullNeZero
import GFullAnalytic
import AnalyticDiffContOnCl

/-!
# §1401 — concrete `hgB` for ζ: `Re(logDeriv g_full (σ+iγ)) ≥ -B` from the factorization data

Threads the ζ-specific suppliers into the two-radii capstone §1378: `hfac` from §1391 (sphere factorization),
`hgnz` from §1382 (inner nonvanishing), `hgdiff` from §1384∘§1383 (`DiffContOnCl` via `R < R'`), `hζsphere` given.
Given the §1343 factorization data `(g₀, hga, hgn, hfac_cod)` on `ball c R'` (`R < R'`), a sphere ζ-bound `Mζ`, and
`hgc_lt`, this produces the Borel–Carathéodory bound for the concrete half-disk factor `g_full = P_rest · g₀`.
Reconciled to §1400's ℤ-exponent `g_full` via §1394, this is §1336's `hgB`.  EXIT0 / std-axioms.
-/

open Complex Metric Classical Filter

theorem hgB_riemannZeta_concrete {c : ℂ} {R' R Mζ : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hball : (1 : ℂ) ∉ ball c R') {g₀ : ℂ → ℂ} (hga : AnalyticOnNhd ℂ g₀ (ball c R'))
    (hgn : ∀ u : ball c R', g₀ u ≠ 0)
    (hfac_cod : riemannZeta =ᶠ[codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u)) • g₀)
    (hζsphere : ∀ s ∈ sphere c R, ‖riemannZeta s‖ ≤ Mζ)
    (hgc_lt : ‖(∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
        (fun u => ¬ u ∈ closedBall c (R / 2)),
        (c - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g₀ c‖
      < Mζ / (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
        (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat)) :
    -(4 * (Real.log (Mζ / (R / 2) ^ (∑ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
        (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat))
        - Real.log ‖(∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
          (fun u => ¬ u ∈ closedBall c (R / 2)),
          (c - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g₀ c‖) / (R / 2))
      ≤ (logDeriv (fun s => (∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
          (fun u => ¬ u ∈ closedBall c (R / 2)),
          (s - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g₀ s) c).re := by
  set m : ℂ → ℕ := fun u => (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat with hm
  set p : ℂ → Prop := fun u => u ∈ closedBall c (R / 2) with hp
  refine re_logDeriv_ge_of_halfBall_factor' (c := c) (R := R) (Mζ := Mζ)
    (P := fun s => ∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter p,
      (s - u) ^ m u)
    (S := (divisor_riemannZeta_ball_support_finite c R').toFinset.filter p) (m := m)
    (g := fun s => (∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
      (fun u => ¬ p u), (s - u) ^ m u) * g₀ s)
    hR rfl ?_ ?_ hζsphere ?_ ?_ hgc_lt
  · -- hS: ∀ρ∈S, ρ∈closedBall c (R/2)
    intro ρ hρ; exact (Finset.mem_filter.mp hρ).2
  · -- hfac: ∀s∈sphere c R, ζ s = P s * g_full s
    intro s hs
    exact zeta_eq_split_on_sphere hRR' hball hga hfac_cod p hs
  · -- hgdiff: DiffContOnCl g_full (ball c R)
    exact analyticOnNhd_diffContOnCl hR hRR' (prod_rest_mul_analyticOnNhd m hga)
  · -- hgnz: ∀z∈ball c (R/2), g_full z ≠ 0
    intro z hz
    refine prod_rest_mul_ne_zero m g₀ ?_ hz (hgn ⟨z, ball_subset_ball (by linarith) hz⟩)
    intro ρ hρ; exact (Finset.mem_filter.mp hρ).2

#print axioms hgB_riemannZeta_concrete
