import OneLeDivisorLFunction
import LFunctionDivisorFinite

/-!
# §1495 — the kept DLVP zero lies in the L divisor support (`χ ≠ 1`)

L-analogue of ζ §1353.  From §1492 the divisor at a genuine zero `ρ = β+iγ` of `LFunction χ` is `≥ 1`, hence
`≠ 0`, so `ρ` lies in `Function.support (divisor (LFunction χ) (ball c R))`, i.e. in the support Finset — the index
set of the region's partial-fractions sum.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem LFunction_zero_mem_divisor_support {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) {c : ℂ} {R : ℝ} {ρ z₀ : ℂ} (hρ : ρ ∈ ball c R) (hLρ : LFunction χ ρ = 0)
    (hz₀ : z₀ ∈ ball c R) (hz₀ne : LFunction χ z₀ ≠ 0) :
    ρ ∈ (divisor_LFunction_ball_support_finite hχ c R).toFinset := by
  rw [Set.Finite.mem_toFinset, Function.mem_support]
  have h1 : 1 ≤ MeromorphicOn.divisor (LFunction χ) (ball c R) ρ :=
    one_le_divisor_LFunction_at_zero hχ hρ hLρ hz₀ hz₀ne
  omega

#print axioms LFunction_zero_mem_divisor_support
