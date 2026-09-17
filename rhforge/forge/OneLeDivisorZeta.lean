import AnalyticAtZeta
import ExplicitFormula
import Mathlib.Analysis.Meromorphic.Divisor
import Mathlib.Analysis.Meromorphic.Order
import Mathlib.Analysis.Analytic.Order
import Mathlib.NumberTheory.LSeries.Dirichlet

/-!
# §1351 — `1 ≤ divisor ζ (ball c R) ρ` at a genuine zero (§1335/§1336's `hm₀`)

On a pole-free disk `B(c,R)` (avoiding `s = 1`) carrying a nonvanishing point `z₀` (`ζ z₀ ≠ 0`), the
kept DLVP zero `ρ` (`ζ ρ = 0`) has divisor `≥ 1`.  Chain: ζ analytic at `ρ` and `z₀` (§1345); order at
the nonvanishing point is `0`, hence `≠ ⊤`; `meromorphicOrderAt_ne_top_of_isPreconnected` propagates
`≠ ⊤` to `ρ` across the connected ball; `AnalyticAt.analyticOrderAt_ne_zero` gives order `≠ 0` from
`ζ ρ = 0`; `AnalyticAt.meromorphicOrderAt_nonneg` gives order `≥ 0`.  A finite integer order that is
`≥ 0` and `≠ 0` is `≥ 1`, so `divisor = order.untop₀ ≥ 1`.  This is §1335/§1336's `hm₀ : 1 ≤ m ρ₀`.
EXIT0 / std-axioms.
-/

open Complex Metric

theorem one_le_divisor_riemannZeta_at_zero {c : ℂ} {R : ℝ} (hball : (1 : ℂ) ∉ ball c R)
    {ρ z₀ : ℂ} (hρ : ρ ∈ ball c R) (hζρ : riemannZeta ρ = 0)
    (hz₀ : z₀ ∈ ball c R) (hz₀ne : riemannZeta z₀ ≠ 0) :
    1 ≤ MeromorphicOn.divisor riemannZeta (ball c R) ρ := by
  have hmero : MeromorphicOn riemannZeta (ball c R) :=
    meromorphicOn_riemannZeta.mono_set (Set.subset_univ _)
  have hρ1 : ρ ≠ 1 := fun h => hball (h ▸ hρ)
  have hz₀1 : z₀ ≠ 1 := fun h => hball (h ▸ hz₀)
  have haρ : AnalyticAt ℂ riemannZeta ρ := analyticAt_riemannZeta hρ1
  have haz₀ : AnalyticAt ℂ riemannZeta z₀ := analyticAt_riemannZeta hz₀1
  -- order at the nonvanishing point z₀ is 0, hence ≠ ⊤
  have hord_z₀ : meromorphicOrderAt riemannZeta z₀ = 0 := by
    rw [haz₀.meromorphicOrderAt_eq, haz₀.analyticOrderAt_eq_zero.mpr hz₀ne]; rfl
  have hz₀top : meromorphicOrderAt riemannZeta z₀ ≠ ⊤ := by
    rw [hord_z₀]; exact WithTop.zero_ne_top
  -- ≠ ⊤ propagates to ρ across the connected ball
  have hρtop : meromorphicOrderAt riemannZeta ρ ≠ ⊤ :=
    hmero.meromorphicOrderAt_ne_top_of_isPreconnected (isPreconnected_ball) hz₀ hρ hz₀top
  -- order at ρ ≠ 0 (since ζ ρ = 0)
  have haor : analyticOrderAt riemannZeta ρ ≠ 0 := haρ.analyticOrderAt_ne_zero.mpr hζρ
  have hρ0 : meromorphicOrderAt riemannZeta ρ ≠ 0 := by
    rw [haρ.meromorphicOrderAt_eq]; simpa using haor
  -- order at ρ ≥ 0 (analytic)
  have hρnn : 0 ≤ meromorphicOrderAt riemannZeta ρ := haρ.meromorphicOrderAt_nonneg
  rw [MeromorphicOn.divisor_apply hmero hρ]
  lift meromorphicOrderAt riemannZeta ρ to ℤ using hρtop with n hn
  have hn0 : n ≠ 0 := by rintro rfl; exact hρ0 (by simpa using hn.symm)
  have hnn : 0 ≤ n := by exact_mod_cast hρnn
  simp only [WithTop.untop₀_coe]; omega

#print axioms one_le_divisor_riemannZeta_at_zero
