import LFunctionFactorization
import LogDerivFromCodiscreteFactor
import LogDerivZpowMulSplit
import ContinuousAtLogDeriv
import FactorAnalyticNeZero
import DivisorLFunctionZeroReGtOne
import SubNeZeroDivisorLFunction
import LFunctionDivisorFinite
import Mathlib.NumberTheory.LSeries.Nonvanishing

/-!
# §1486 — half-disk logDeriv decomposition of `LFunction χ` with **external** factorization `g`

L-analogue of the ζ coordination brick §1400.  Takes the §1476 factorization data `(g, hga, hgn, hfac)` as
hypotheses so the final assembly obtains §1476 **once** and shares the same `g` (hence the same `g_full`) across
both the half-disk decomposition (this brick, feeding `hdecomp`) and the two-radii BC bound (§1481/§1482, feeding
`hgB`) — resolving the g-coordination.  For a test point `z` (`Re z > 1`) and predicate `p`:

  `logDeriv (LFunction χ) z = ∑_{ρ ∈ filter p} divisor ρ/(z-ρ) + logDeriv g_full z`,
  `g_full = (∏_{ρ ∈ filter ¬p} (·-ρ)^{divisor ρ}) · g`.

Proof = §1485's fact-discharge + generic §1342 + Finset split + generic §1332 (regroup annulus into `g_full`).
EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Filter Classical

theorem logDeriv_LFunction_halfdisk_decomp_of_factor {N : ℕ} [NeZero N]
    {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1) {c : ℂ} {R : ℝ} {z : ℂ}
    (hz : z ∈ ball c R) (hzre : 1 < z.re) {g : ℂ → ℂ}
    (hga : AnalyticOnNhd ℂ g (ball c R)) (hgn : ∀ u : ball c R, g u ≠ 0)
    (hfac : LFunction χ =ᶠ[codiscreteWithin (ball c R)]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R) u)) • g)
    (p : ℂ → Prop) :
    logDeriv (LFunction χ) z =
      (∑ ρ ∈ (divisor_LFunction_ball_support_finite hχ c R).toFinset.filter p,
        (MeromorphicOn.divisor (LFunction χ) (ball c R) ρ : ℂ) / (z - ρ))
      + logDeriv (fun s => (∏ ρ ∈ (divisor_LFunction_ball_support_finite hχ c R).toFinset.filter
          (fun u => ¬ p u), (s - ρ) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R) ρ)) * g s) z := by
  set T := (divisor_LFunction_ball_support_finite hχ c R).toFinset with hT
  set D := MeromorphicOn.divisor (LFunction χ) (ball c R) with hD
  have hana : AnalyticOnNhd ℂ (LFunction χ) Set.univ :=
    (differentiable_LFunction hχ).differentiableOn.analyticOnNhd isOpen_univ
  have hLz : LFunction χ z ≠ 0 := LFunction_ne_zero_of_one_le_re χ (Or.inl hχ) hzre.le
  have hζc : ContinuousAt (logDeriv (LFunction χ)) z :=
    continuousAt_logDeriv (hana z (Set.mem_univ z)) hLz
  have hdz : D z = 0 := divisor_LFunction_eq_zero_of_one_lt_re hχ hz hzre
  have hga_z : AnalyticAt ℂ g z := hga z hz
  have hgz : g z ≠ 0 := hgn ⟨z, hz⟩
  obtain ⟨hFa, hFne⟩ := analyticAt_and_ne_zero_factor hdz hga_z hgz
  have hFc : ContinuousAt (logDeriv ((∏ᶠ u, (· - u) ^ (D u)) • g)) z :=
    continuousAt_logDeriv hFa hFne
  have hsub : ∀ ρ ∈ T, z - ρ ≠ 0 := by
    intro ρ hρ
    rw [hT, Set.Finite.mem_toFinset, Function.mem_support] at hρ
    exact sub_ne_zero_of_divisor_LFunction_ne_zero hχ (D.supportWithinDomain hρ) hρ hzre
  have hfull := logDeriv_eq_sum_of_codiscrete_factor isOpen_ball hz
    (divisor_LFunction_ball_support_finite hχ c R) hfac hζc hFc hsub hgz hga_z.differentiableAt
  have hsub' : ∀ ρ ∈ T.filter (fun u => ¬ p u), z - ρ ≠ 0 :=
    fun ρ hρ => hsub ρ (Finset.mem_of_mem_filter ρ hρ)
  have hsplit := logDeriv_prod_zpow_mul_split (S := T.filter (fun u => ¬ p u)) (m := D)
    (h := g) (z := z) hsub' hgz hga_z.differentiableAt
  rw [hfull, ← Finset.sum_filter_add_sum_filter_not T p, hsplit]
  ring

#print axioms logDeriv_LFunction_halfdisk_decomp_of_factor
