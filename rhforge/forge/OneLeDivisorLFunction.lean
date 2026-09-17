import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.Analysis.Meromorphic.Divisor
import Mathlib.Analysis.Meromorphic.Order
import Mathlib.Analysis.Analytic.Order
import Mathlib.Analysis.Complex.CauchyIntegral

/-!
# §1492 — `1 ≤ divisor (LFunction χ) (ball c R) ρ` at a genuine zero (`χ ≠ 1`; the `hm₀` input)

L-analogue of ζ §1351 (cleaner: `LFunction χ` entire, no pole exception).  On a ball carrying a nonvanishing
point `z₀`, the kept DLVP zero `ρ` (`LFunction χ ρ = 0`) has divisor `≥ 1`.  Chain: `LFunction χ` analytic at `ρ`
and `z₀`; order at `z₀` is `0` (nonvanishing) hence `≠ ⊤`; `meromorphicOrderAt_ne_top_of_isPreconnected`
propagates `≠ ⊤` to `ρ`; `analyticOrderAt_ne_zero` gives order `≠ 0` from `LFunction χ ρ = 0`; nonneg from analytic.
A finite integer order `≥ 0` and `≠ 0` is `≥ 1`, so `divisor ≥ 1`.  This is the `hm₀ : 1 ≤ m ρ₀` the keep-one-zero
step §1335 and the region composition §1336 consume.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem one_le_divisor_LFunction_at_zero {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c : ℂ} {R : ℝ} {ρ z₀ : ℂ} (hρ : ρ ∈ ball c R) (hLρ : LFunction χ ρ = 0)
    (hz₀ : z₀ ∈ ball c R) (hz₀ne : LFunction χ z₀ ≠ 0) :
    1 ≤ MeromorphicOn.divisor (LFunction χ) (ball c R) ρ := by
  have hana : AnalyticOnNhd ℂ (LFunction χ) Set.univ :=
    (differentiable_LFunction hχ).differentiableOn.analyticOnNhd isOpen_univ
  have hmero : MeromorphicOn (LFunction χ) (ball c R) :=
    hana.meromorphicOn.mono_set (Set.subset_univ _)
  have haρ : AnalyticAt ℂ (LFunction χ) ρ := hana ρ (Set.mem_univ ρ)
  have haz₀ : AnalyticAt ℂ (LFunction χ) z₀ := hana z₀ (Set.mem_univ z₀)
  have hord_z₀ : meromorphicOrderAt (LFunction χ) z₀ = 0 := by
    rw [haz₀.meromorphicOrderAt_eq, haz₀.analyticOrderAt_eq_zero.mpr hz₀ne]; rfl
  have hz₀top : meromorphicOrderAt (LFunction χ) z₀ ≠ ⊤ := by rw [hord_z₀]; exact WithTop.zero_ne_top
  have hρtop : meromorphicOrderAt (LFunction χ) ρ ≠ ⊤ :=
    hmero.meromorphicOrderAt_ne_top_of_isPreconnected isPreconnected_ball hz₀ hρ hz₀top
  have haor : analyticOrderAt (LFunction χ) ρ ≠ 0 := haρ.analyticOrderAt_ne_zero.mpr hLρ
  have hρ0 : meromorphicOrderAt (LFunction χ) ρ ≠ 0 := by
    rw [haρ.meromorphicOrderAt_eq]; simpa using haor
  have hρnn : 0 ≤ meromorphicOrderAt (LFunction χ) ρ := haρ.meromorphicOrderAt_nonneg
  rw [MeromorphicOn.divisor_apply hmero hρ]
  lift meromorphicOrderAt (LFunction χ) ρ to ℤ using hρtop with n hn
  have hn0 : n ≠ 0 := by rintro rfl; exact hρ0 (by simpa using hn.symm)
  have hnn : 0 ≤ n := by exact_mod_cast hρnn
  simp only [WithTop.untop₀_coe]; omega

#print axioms one_le_divisor_LFunction_at_zero
