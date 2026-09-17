import ZetaDivisorNonneg

/-!
# §1387 — `0 ≤ divisor ζ (ball c R) u` for **all** `u` (pole-free ball)

The `∀ u` divisor-nonnegativity §1386 needs: inside the pole-free ball it is §1350 (`divisor ≥ 0`); outside the
ball `u` is not in the (within-domain) support, so `divisor u = 0`.  Either way `0 ≤ divisor u`.  This discharges
`smul_split_eval_nat`'s `hdnn` hypothesis when specialising the split-evaluation to ζ.  EXIT0 / std-axioms.
-/

open Metric

theorem divisor_riemannZeta_nonneg_all {c : ℂ} {R : ℝ} (hball : (1 : ℂ) ∉ ball c R) (u : ℂ) :
    0 ≤ MeromorphicOn.divisor riemannZeta (ball c R) u := by
  by_cases hu : u ∈ ball c R
  · exact divisor_riemannZeta_nonneg hball hu
  · have hns : u ∉ Function.support (MeromorphicOn.divisor riemannZeta (ball c R)) :=
      fun h => hu ((MeromorphicOn.divisor riemannZeta (ball c R)).supportWithinDomain h)
    exact le_of_eq (Function.notMem_support.mp hns).symm

#print axioms divisor_riemannZeta_nonneg_all
