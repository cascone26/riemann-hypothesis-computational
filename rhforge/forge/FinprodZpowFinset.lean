import Mathlib.Algebra.BigOperators.Finprod
import Mathlib.Analysis.SpecialFunctions.Complex.Analytic

/-!
# §1337 — `∏ᶠ u, (x-u)^(d u) = ∏ u ∈ support, (x-u)^(d u)` (finprod → finset conversion)

The reusable conversion turning Mathlib's factorized-rational finprod `∏ᶠ u, (·-u)^(d u)` (as `extract_zeros_poles`
produces, `d` the ζ-divisor) into a finite product over the divisor support.  Factors with `d u = 0` are
`(x-u)^0 = 1`, so `mulSupport (fun u => (x-u)^(d u)) ⊆ support d`, and `finprod_eq_prod_of_mulSupport_subset` gives
the finite product.  This is the bridge from the extract_zeros_poles output to the finset form §1333/§1336 consume.
EXIT0 / std-axioms.
-/

open Function

theorem finprod_zpow_sub_eq_finset_prod {d : ℂ → ℤ} {x : ℂ} (hd : (Function.support d).Finite) :
    (∏ᶠ u, (x - u) ^ (d u)) = ∏ u ∈ hd.toFinset, (x - u) ^ (d u) := by
  apply finprod_eq_prod_of_mulSupport_subset
  intro u hu
  rw [Set.Finite.coe_toFinset, Function.mem_support]
  intro hdu
  exact hu (by simp [hdu])

#print axioms finprod_zpow_sub_eq_finset_prod
