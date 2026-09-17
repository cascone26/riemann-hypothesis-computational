import DivisorLFunctionZeroReGtOne

/-!
# §1484 — Test point avoids the L divisor support: `z - ρ ≠ 0` for support points `ρ`

L-analogue of ζ §1348.  Any `ρ ∈ U` with `divisor (LFunction χ) U ρ ≠ 0` (a zero of `LFunction χ`) has
`Re ρ ≤ 1` (contrapositive of §1483), so for the DLVP test point `z` (`Re z > 1`), `Re z > 1 ≥ Re ρ` forces
`z ≠ ρ`, i.e. `z - ρ ≠ 0`.  This is the `hsub` hypothesis of the generic decomposition §1342 for `LFunction χ`.
EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem sub_ne_zero_of_divisor_LFunction_ne_zero {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) {U : Set ℂ} {z ρ : ℂ} (hρU : ρ ∈ U)
    (hρ : MeromorphicOn.divisor (LFunction χ) U ρ ≠ 0) (hzre : 1 < z.re) : z - ρ ≠ 0 := by
  have hρre : ρ.re ≤ 1 := by
    by_contra h
    exact hρ (divisor_LFunction_eq_zero_of_one_lt_re hχ hρU (not_le.mp h))
  intro hzρ
  rw [sub_eq_zero] at hzρ
  rw [hzρ] at hzre
  linarith

#print axioms sub_ne_zero_of_divisor_LFunction_ne_zero
