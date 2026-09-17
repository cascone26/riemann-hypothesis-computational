import ZetaFactorization
import LogDerivFromCodiscreteFactor
import ContinuousAtLogDeriv
import AnalyticAtZeta
import FactorAnalyticNeZero
import ZetaDivisorZeroReGtOne
import SubNeZeroDivisor
import ZetaDivisorFinite
import Mathlib.NumberTheory.LSeries.Dirichlet

/-!
# §1349 — `logDeriv ζ` decomposition at the test point (the FINAL COMPOSE)

Threading the ζ factorization §1343 through the decomposition assembly §1342, supplying every hypothesis from the
now-proven test-point facts (§1344 continuity, §1345 ζ analytic, §1346 factor analytic≠0, §1347 divisor=0 at Re>1,
§1348 z-ρ≠0).  For a disk `B(c,R)` with a nonvanishing point `z₀` and a test point `z` with `Re z > 1`:

  `∃ g, g analytic nonvanishing on B(c,R) ∧ logDeriv ζ z = ∑_{ρ ∈ support} divisor ρ/(z-ρ) + logDeriv g z`.

This is the concrete `hdecomp` of ζ obtained entirely from Mathlib + the verified converters — the exact input
§1335/§1336 turn into the quantitative zero-free region.  EXIT0 / std-axioms.
-/

open Complex Metric Filter

theorem logDeriv_riemannZeta_decomp_at {c : ℂ} {R : ℝ} {z z₀ : ℂ}
    (hz : z ∈ ball c R) (hzre : 1 < z.re)
    (hz₀ : z₀ ∈ ball c R) (hz₀1 : z₀ ≠ 1) (hz₀ne : riemannZeta z₀ ≠ 0) :
    ∃ g : ℂ → ℂ, AnalyticOnNhd ℂ g (ball c R) ∧ (∀ u : ball c R, g u ≠ 0) ∧
      logDeriv riemannZeta z =
        (∑ ρ ∈ (divisor_riemannZeta_ball_support_finite c R).toFinset,
          (MeromorphicOn.divisor riemannZeta (ball c R) ρ : ℂ) / (z - ρ)) + logDeriv g z := by
  obtain ⟨g, hga, hgn, hfac⟩ := exists_factorization_riemannZeta hz₀ hz₀1 hz₀ne
  have hz1 : z ≠ 1 := fun h => by rw [h] at hzre; simp at hzre
  have hζc : ContinuousAt (logDeriv riemannZeta) z :=
    continuousAt_logDeriv (analyticAt_riemannZeta hz1) (riemannZeta_ne_zero_of_one_lt_re hzre)
  have hdz : MeromorphicOn.divisor riemannZeta (ball c R) z = 0 :=
    divisor_riemannZeta_eq_zero_of_one_lt_re hz hzre
  have hga_z : AnalyticAt ℂ g z := hga z hz
  have hgz : g z ≠ 0 := hgn ⟨z, hz⟩
  obtain ⟨hFa, hFne⟩ := analyticAt_and_ne_zero_factor hdz hga_z hgz
  have hFc : ContinuousAt
      (logDeriv ((∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R) u)) • g)) z :=
    continuousAt_logDeriv hFa hFne
  have hsub : ∀ ρ ∈ (divisor_riemannZeta_ball_support_finite c R).toFinset, z - ρ ≠ 0 := by
    intro ρ hρ
    rw [Set.Finite.mem_toFinset, Function.mem_support] at hρ
    have hρU : ρ ∈ ball c R := (MeromorphicOn.divisor riemannZeta (ball c R)).supportWithinDomain hρ
    exact sub_ne_zero_of_divisor_ne_zero hρU hρ hzre
  refine ⟨g, hga, hgn, ?_⟩
  exact logDeriv_eq_sum_of_codiscrete_factor isOpen_ball hz
    (divisor_riemannZeta_ball_support_finite c R) hfac hζc hFc hsub hgz hga_z.differentiableAt

#print axioms logDeriv_riemannZeta_decomp_at
