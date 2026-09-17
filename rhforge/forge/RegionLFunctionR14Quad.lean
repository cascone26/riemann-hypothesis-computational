import RegionLFunctionOrFarQuad

set_option maxHeartbeats 1000000

/-!
# §1560 — QUADRATIC L region with geometry fixed (`R = 1/4`, `R' = |γ|`)

Quadratic analogue of §1545.  Specializes §1559 to `R = 1/4`, `R' = |γ|`; `hfar` via §1434 (`K ≥ 32/7`), the
radius conditions reduce to `1 < |γ|` (using `|2γ| = 2|γ|`).  Region `β ≤ 1-1/(14K)` for `LFunction χ` (χ≠1,
χ²=1) from just: the zero, `1 < |γ|`, `K ≥ 32/7`, `A,B,C ≥ 0`, the two log bounds (at R=1/4), `hA_pole`.
EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem region_LFunction_R14_quad {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) (hχ2eq : χ ^ 2 = 1) {γ β A B C : ℝ}
    (hγ : 1 < |γ|) (hKbig : 32 / 7 ≤ 3 * A + 4 * B + C)
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hLρ₀ : LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hBexp : 4 * Real.log ((2 + (N : ℝ)
        + ((1 + 1 / (2 * (3 * A + 4 * B + C))) + |γ| + 1 / 4) * (N : ℝ)
          / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1 / 4))
        * (3 + 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1))) / (1 / 4 / 2) ≤ B)
    (hCexp : 4 * Real.log ((3 + 1 / (|2 * γ| - 1 / 4)
        + ((1 + 1 / (2 * (3 * A + 4 * B + C))) + |2 * γ| + 1 / 4)
          / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1 / 4))
        * (3 + 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1))) / (1 / 4 / 2)
      ≤ C - ∑ p ∈ N.primeFactors, Real.log p)
    (hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ)) : ℂ))‖ ≤ A) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hK : 0 < 3 * A + 4 * B + C := by linarith
  have h2γ : |2 * γ| = 2 * |γ| := by rw [abs_mul, abs_two]
  have hfar : (1 + 1 / (2 * (3 * A + 4 * B + C))) - (1 / 4 : ℝ) / 2
      ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) :=
    hfar_holds hK (by nlinarith [hKbig])
  exact region_LFunction_or_far_quad (R := 1 / 4) (R' := |γ|) hχ hχ2eq (by norm_num) (by norm_num)
    (by linarith) (by rw [h2γ]; linarith) (by rw [h2γ]; linarith [abs_nonneg γ])
    hA hB hC hK hfar hLρ₀ hBexp hCexp hA_pole

#print axioms region_LFunction_R14_quad
