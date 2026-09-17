import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.Analysis.Meromorphic.Divisor
import Mathlib.Analysis.Meromorphic.Order
import Mathlib.Analysis.Complex.CauchyIntegral

/-!
# §1479 — `0 ≤ divisor (LFunction χ) (ball c R) u` for **all** `u` (`χ ≠ 1`)

The L-analogue of ζ §1350+§1387, but cleaner: `LFunction χ` is *entire* (`differentiable_LFunction`), so it is
analytic at every point with no pole exception.  Inside the ball, `AnalyticAt.meromorphicOrderAt_nonneg` gives
`order ≥ 0`, hence `divisor = order.untop₀ ≥ 0`; outside the ball `u ∉ support`, so `divisor u = 0`.  This
discharges the `hdnn` (divisor ≥ 0) hypothesis of the generic split-evaluation `smul_split_eval_nat` for
`LFunction χ`, the bridge toward the `P_half · g_full` sphere split.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem divisor_LFunction_nonneg_all {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c : ℂ} {R : ℝ} (u : ℂ) :
    0 ≤ MeromorphicOn.divisor (LFunction χ) (ball c R) u := by
  have hana : AnalyticOnNhd ℂ (LFunction χ) Set.univ :=
    (differentiable_LFunction hχ).differentiableOn.analyticOnNhd isOpen_univ
  have hmero : MeromorphicOn (LFunction χ) (ball c R) :=
    hana.meromorphicOn.mono_set (Set.subset_univ _)
  by_cases hu : u ∈ ball c R
  · have hord : 0 ≤ meromorphicOrderAt (LFunction χ) u :=
      (hana u (Set.mem_univ u)).meromorphicOrderAt_nonneg
    rw [MeromorphicOn.divisor_apply hmero hu]
    rcases eq_or_ne (meromorphicOrderAt (LFunction χ) u) ⊤ with h | h
    · rw [h]; simp
    · lift meromorphicOrderAt (LFunction χ) u to ℤ using h with n
      simpa using hord
  · have hns : u ∉ Function.support (MeromorphicOn.divisor (LFunction χ) (ball c R)) :=
      fun h => hu ((MeromorphicOn.divisor (LFunction χ) (ball c R)).supportWithinDomain h)
    exact le_of_eq (Function.notMem_support.mp hns).symm

#print axioms divisor_LFunction_nonneg_all
