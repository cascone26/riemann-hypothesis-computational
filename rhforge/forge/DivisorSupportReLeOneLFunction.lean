import DivisorLFunctionZeroReGtOne
import LFunctionDivisorFinite

/-!
# §1489 — L divisor-support points satisfy `Re ρ ≤ 1` (`χ ≠ 1`)

L-analogue of ζ §1354.  Every point `ρ` of the `LFunction χ` divisor-support on `ball c R` (a zero of
`LFunction χ`) has `Re ρ ≤ 1`: contrapositive of §1483 (`divisor = 0` whenever `Re > 1`).  Feeds the filtered
`hzre` (`Re ρ ≤ Re z`) the drop-zeros growth bound §1399 needs.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem divisor_support_re_le_one_LFunction {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) {c : ℂ} {R : ℝ} {ρ : ℂ}
    (hρmem : ρ ∈ (divisor_LFunction_ball_support_finite hχ c R).toFinset) : ρ.re ≤ 1 := by
  rw [Set.Finite.mem_toFinset, Function.mem_support] at hρmem
  by_contra h
  have hρU : ρ ∈ ball c R :=
    (MeromorphicOn.divisor (LFunction χ) (ball c R)).supportWithinDomain hρmem
  exact hρmem (divisor_LFunction_eq_zero_of_one_lt_re hχ hρU (not_le.mp h))

#print axioms divisor_support_re_le_one_LFunction
