import RegionLFunctionR14
import HBexpRhsLeLFunction

set_option maxHeartbeats 1000000

/-!
# §1547 — quantitative L(·,χ) region from CLEAN linear-in-γ log bounds

L-analogue of ζ §1457.  Discharges §1545's giant `hBexp`/`hCexp` from the clean bounds
`B ≥ 32·log(N·(2|γ|+6)·(3+2K))` and `C ≥ 32·log(N·(2|2γ|+6)·(3+2K))` (via §1546 + `le_trans`).  The quantitative
DLVP region `β ≤ 1-1/(14K)` for `LFunction χ` (χ≠1, χ²≠1) now follows from: the zero, `1<|γ|`, `K ≥ 32/7`,
`A,B,C ≥ 0`, `hA_pole`, and these two linear-in-`|γ|` log bounds — the shape the §1445 threshold closes.  `N ≥ 1`
from `NeZero N`.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem region_LFunction_R14_clean {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) (hχ2 : χ ^ 2 ≠ 1) {γ β A B C : ℝ}
    (hγ : 1 < |γ|) (hKbig : 32 / 7 ≤ 3 * A + 4 * B + C)
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hLρ₀ : LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hBclean : 32 * Real.log ((N : ℝ) * (2 * |γ| + 6) * (3 + 2 * (3 * A + 4 * B + C))) ≤ B)
    (hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ)) : ℂ))‖ ≤ A)
    (hCclean : 32 * Real.log ((N : ℝ) * (2 * |2 * γ| + 6) * (3 + 2 * (3 * A + 4 * B + C))) ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hKhalf : (1 : ℝ) / 2 ≤ 3 * A + 4 * B + C := by linarith
  have hγ2 : (1 : ℝ) < |2 * γ| := by rw [abs_mul]; simp only [abs_two]; linarith [abs_nonneg γ]
  have hN1 : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast Nat.one_le_iff_ne_zero.mpr (NeZero.ne N)
  exact region_LFunction_R14 hχ hχ2 hγ hKbig hA hB hC hLρ₀
    (le_trans (hBexp_rhs_le_LFunction hKhalf hγ hN1) hBclean)
    (le_trans (hBexp_rhs_le_LFunction hKhalf hγ2 hN1) hCclean) hA_pole

#print axioms region_LFunction_R14_clean
