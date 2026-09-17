import RegionLFunctionOrFar

set_option maxHeartbeats 1000000

/-!
# §1545 — L region with geometry fixed (`R = 1/4`, `R' = |γ|`)

L-analogue of ζ §1454.  Specializes §1544 to `R = 1/4`, `R' = |γ|`.  `hfar` is discharged by §1434 (`8 ≤ 7K/4`,
i.e. `K ≥ 32/7`), and `hRR'` reduces to `1/4 < |γ|`.  The region `β ≤ 1-1/(14K)` for `LFunction χ` (χ≠1, χ²≠1) now
follows from just: the zero, `1 < |γ|`, `K ≥ 32/7`, `A,B,C ≥ 0`, the two explicit log bounds `hBexp`/`hCexp`
(at `R=1/4`), and `hA_pole`.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem region_LFunction_R14 {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) (hχ2 : χ ^ 2 ≠ 1) {γ β A B C : ℝ}
    (hγ : 1 < |γ|) (hKbig : 32 / 7 ≤ 3 * A + 4 * B + C)
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hLρ₀ : LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hBexp : 4 * Real.log ((2 + (N : ℝ)
        + ((1 + 1 / (2 * (3 * A + 4 * B + C))) + |γ| + 1 / 4) * (N : ℝ)
          / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1 / 4))
        * (3 + 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1))) / (1 / 4 / 2) ≤ B)
    (hCexp : 4 * Real.log ((2 + (N : ℝ)
        + ((1 + 1 / (2 * (3 * A + 4 * B + C))) + |2 * γ| + 1 / 4) * (N : ℝ)
          / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1 / 4))
        * (3 + 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1))) / (1 / 4 / 2) ≤ C)
    (hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ)) : ℂ))‖ ≤ A) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hK : 0 < 3 * A + 4 * B + C := by linarith
  have hfar : (1 + 1 / (2 * (3 * A + 4 * B + C))) - (1 / 4 : ℝ) / 2
      ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) :=
    hfar_holds hK (by nlinarith [hKbig])
  exact region_LFunction_or_far (R := 1 / 4) (R' := |γ|) hχ hχ2 (by norm_num) (by norm_num)
    (by linarith) hA hB hC hK hfar hLρ₀ hBexp hCexp hA_pole

#print axioms region_LFunction_R14
