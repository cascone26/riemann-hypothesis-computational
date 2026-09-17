import ZetaDivisorZeroReGtOne

/-!
# §1348 — Test point avoids the divisor support: `z - ρ ≠ 0` for support points `ρ`

Any `ρ ∈ U` with `divisor ζ U ρ ≠ 0` (a zero or pole of ζ) has `Re ρ ≤ 1` (contrapositive of §1347), so for the DLVP
test point `z` (`Re z > 1`), `Re z > 1 ≥ Re ρ` forces `z ≠ ρ`, i.e. `z - ρ ≠ 0`.  This is §1342's `hsub` hypothesis
(`∀ ρ ∈ divisor-support, z - ρ ≠ 0`).  EXIT0 / std-axioms.
-/

open Complex

theorem sub_ne_zero_of_divisor_ne_zero {U : Set ℂ} {z ρ : ℂ} (hρU : ρ ∈ U)
    (hρ : MeromorphicOn.divisor riemannZeta U ρ ≠ 0) (hzre : 1 < z.re) : z - ρ ≠ 0 := by
  have hρre : ρ.re ≤ 1 := by
    by_contra h
    exact hρ (divisor_riemannZeta_eq_zero_of_one_lt_re hρU (not_le.mp h))
  intro hzρ
  rw [sub_eq_zero] at hzρ
  rw [hzρ] at hzre
  linarith

#print axioms sub_ne_zero_of_divisor_ne_zero
