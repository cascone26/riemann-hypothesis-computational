import RegionLFunctionR14Quad
import HBexpRhsLeLFunction
import HBexpRhsLe

set_option maxHeartbeats 1000000

/-!
# §1561 — QUADRATIC L region from CLEAN linear-in-γ log bounds

Quadratic analogue of §1547.  Discharges §1560's `hBexp` (χ side) via §1546 (L-form → `32·log(N(2|γ|+6)(3+2K))`)
and `hCexp` (χ² side, ζ-form) via §1456 `hBexp_rhs_le` (already proven in the ζ tower → `32·log((3|2γ|+10)(3+2K))`).
The region follows from the clean bounds `hBclean` (χ) and `hCclean` (χ², ζ-form `≤ C - ∑_{p|N} log p`).  `N ≥ 1`
from `NeZero`.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem region_LFunction_R14_clean_quad {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) (hχ2eq : χ ^ 2 = 1) {γ β A B C : ℝ}
    (hγ : 1 < |γ|) (hKbig : 32 / 7 ≤ 3 * A + 4 * B + C)
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hLρ₀ : LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hBclean : 32 * Real.log ((N : ℝ) * (2 * |γ| + 6) * (3 + 2 * (3 * A + 4 * B + C))) ≤ B)
    (hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ)) : ℂ))‖ ≤ A)
    (hCclean : 32 * Real.log ((3 * |2 * γ| + 10) * (3 + 2 * (3 * A + 4 * B + C)))
      ≤ C - ∑ p ∈ N.primeFactors, Real.log p) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hKhalf : (1 : ℝ) / 2 ≤ 3 * A + 4 * B + C := by linarith
  have hγ2 : (1 : ℝ) < |2 * γ| := by rw [abs_mul]; simp only [abs_two]; linarith [abs_nonneg γ]
  have hN1 : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast Nat.one_le_iff_ne_zero.mpr (NeZero.ne N)
  exact region_LFunction_R14_quad hχ hχ2eq hγ hKbig hA hB hC hLρ₀
    (le_trans (hBexp_rhs_le_LFunction hKhalf hγ hN1) hBclean)
    (le_trans (hBexp_rhs_le hKhalf hγ2) hCclean) hA_pole

#print axioms region_LFunction_R14_clean_quad
