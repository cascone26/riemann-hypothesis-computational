import ReLogDerivGeHalfBallGen
import LFunctionEqSplitSphere
import GFullNeZero
import GFullAnalytic
import AnalyticDiffContOnCl

/-!
# §1482 — concrete `hgB` for `LFunction χ`: `Re(logDeriv g_full (σ+iγ)) ≥ -B` from factorization data

L-analogue of the ζ capstone §1401, now built on the GENERIC Borel–Carathéodory factor bound §1471
(`re_logDeriv_ge_of_halfBall_factor_gen`).  Threads the L suppliers into it: `hfac` from §1481 (L sphere-split),
`hgnz` from §1382 (inner nonvanishing, generic), `hgdiff` from §1384∘§1383 (`DiffContOnCl` via `R < R'`, generic),
`hFsphere` given.  Cleaner than ζ — `LFunction χ` is entire, so no `1 ∉ ball` pole hypothesis.  Given the §1476
factorization data `(g₀, hga, hgn, hfac_cod)` on `ball c R'` (`R < R'`), a sphere L-bound `Mζ`, and `hgc_lt`, this
produces the BC bound for the concrete half-disk factor `g_full = P_rest · g₀` — the L analogue of §1336's `hgB`,
discharging the `Q` term of the L zero-free region §1470.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Classical Filter

theorem hgB_LFunction_concrete {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c : ℂ} {R' R Mζ : ℝ} (hR : 0 < R) (hRR' : R < R')
    {g₀ : ℂ → ℂ} (hga : AnalyticOnNhd ℂ g₀ (ball c R'))
    (hgn : ∀ u : ball c R', g₀ u ≠ 0)
    (hfac_cod : LFunction χ =ᶠ[codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u)) • g₀)
    (hLsphere : ∀ s ∈ sphere c R, ‖LFunction χ s‖ ≤ Mζ)
    (hgc_lt : ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
        (fun u => ¬ u ∈ closedBall c (R / 2)),
        (c - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g₀ c‖
      < Mζ / (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
        (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat)) :
    -(4 * (Real.log (Mζ / (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
        (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat))
        - Real.log ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
          (fun u => ¬ u ∈ closedBall c (R / 2)),
          (c - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g₀ c‖) / (R / 2))
      ≤ (logDeriv (fun s => (∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
          (fun u => ¬ u ∈ closedBall c (R / 2)),
          (s - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g₀ s) c).re := by
  set m : ℂ → ℕ := fun u => (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat with hm
  set p : ℂ → Prop := fun u => u ∈ closedBall c (R / 2) with hp
  refine re_logDeriv_ge_of_halfBall_factor_gen (F := LFunction χ) (c := c) (R := R) (Mζ := Mζ)
    (P := fun s => ∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter p,
      (s - u) ^ m u)
    (S := (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter p) (m := m)
    (g := fun s => (∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
      (fun u => ¬ p u), (s - u) ^ m u) * g₀ s)
    hR rfl ?_ ?_ hLsphere ?_ ?_ hgc_lt
  · intro ρ hρ; exact (Finset.mem_filter.mp hρ).2
  · intro s hs
    exact LFunction_eq_split_on_sphere hχ hRR' hga hfac_cod p hs
  · exact analyticOnNhd_diffContOnCl hR hRR' (prod_rest_mul_analyticOnNhd m hga)
  · intro z hz
    refine prod_rest_mul_ne_zero m g₀ ?_ hz (hgn ⟨z, ball_subset_ball (by linarith) hz⟩)
    intro ρ hρ; exact (Finset.mem_filter.mp hρ).2

#print axioms hgB_LFunction_concrete
