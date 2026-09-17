import LFunctionFactorization
import LogDerivFromCodiscreteFactor
import ContinuousAtLogDeriv
import FactorAnalyticNeZero
import DivisorLFunctionZeroReGtOne
import SubNeZeroDivisorLFunction
import LFunctionDivisorFinite
import Mathlib.NumberTheory.LSeries.Nonvanishing

/-!
# §1485 — `logDeriv (LFunction χ)` decomposition at the test point (the FINAL COMPOSE)

L-analogue of the ζ capstone §1349.  Threading the L factorization §1476 through the generic decomposition
assembly §1342, supplying every hypothesis from the now-proven L test-point facts (generic §1344 continuity,
§1473 analytic, generic §1346 factor analytic≠0, §1483 divisor=0 at Re>1, §1484 z−ρ≠0).  For a disk `B(c,R)` with
a nonvanishing point `z₀` and a test point `z` with `Re z > 1`:

  `∃ g, g analytic nonvanishing on B(c,R) ∧ logDeriv (LFunction χ) z = ∑_{ρ} divisor ρ/(z-ρ) + logDeriv g z`.

This is the concrete `hdecomp` of `LFunction χ` — the exact `-L'/L` decomposition the DLVP region argument turns
into the quantitative zero-free region, now with the zero-sum `∑ divisor ρ/(z-ρ)` explicit.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Filter

theorem logDeriv_LFunction_decomp_at {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c : ℂ} {R : ℝ} {z z₀ : ℂ} (hz : z ∈ ball c R) (hzre : 1 < z.re)
    (hz₀ : z₀ ∈ ball c R) (hz₀ne : LFunction χ z₀ ≠ 0) :
    ∃ g : ℂ → ℂ, AnalyticOnNhd ℂ g (ball c R) ∧ (∀ u : ball c R, g u ≠ 0) ∧
      logDeriv (LFunction χ) z =
        (∑ ρ ∈ (divisor_LFunction_ball_support_finite hχ c R).toFinset,
          (MeromorphicOn.divisor (LFunction χ) (ball c R) ρ : ℂ) / (z - ρ)) + logDeriv g z := by
  have hana : AnalyticOnNhd ℂ (LFunction χ) Set.univ :=
    (differentiable_LFunction hχ).differentiableOn.analyticOnNhd isOpen_univ
  obtain ⟨g, hga, hgn, hfac⟩ := exists_factorization_LFunction hχ hz₀ hz₀ne
  have hLz : LFunction χ z ≠ 0 := LFunction_ne_zero_of_one_le_re χ (Or.inl hχ) hzre.le
  have hζc : ContinuousAt (logDeriv (LFunction χ)) z :=
    continuousAt_logDeriv (hana z (Set.mem_univ z)) hLz
  have hdz : MeromorphicOn.divisor (LFunction χ) (ball c R) z = 0 :=
    divisor_LFunction_eq_zero_of_one_lt_re hχ hz hzre
  have hga_z : AnalyticAt ℂ g z := hga z hz
  have hgz : g z ≠ 0 := hgn ⟨z, hz⟩
  obtain ⟨hFa, hFne⟩ := analyticAt_and_ne_zero_factor hdz hga_z hgz
  have hFc : ContinuousAt
      (logDeriv ((∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R) u)) • g)) z :=
    continuousAt_logDeriv hFa hFne
  have hsub : ∀ ρ ∈ (divisor_LFunction_ball_support_finite hχ c R).toFinset, z - ρ ≠ 0 := by
    intro ρ hρ
    rw [Set.Finite.mem_toFinset, Function.mem_support] at hρ
    have hρU : ρ ∈ ball c R :=
      (MeromorphicOn.divisor (LFunction χ) (ball c R)).supportWithinDomain hρ
    exact sub_ne_zero_of_divisor_LFunction_ne_zero hχ hρU hρ hzre
  refine ⟨g, hga, hgn, ?_⟩
  exact logDeriv_eq_sum_of_codiscrete_factor isOpen_ball hz
    (divisor_LFunction_ball_support_finite hχ c R) hfac hζc hFc hsub hgz hga_z.differentiableAt

#print axioms logDeriv_LFunction_decomp_at
