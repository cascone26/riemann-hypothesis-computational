import NormLFunctionGeExplicit
import Mathlib.NumberTheory.LSeries.Nonvanishing

set_option maxHeartbeats 1000000

/-!
# §1539 — discharge the BC scalar `hBexplicit` for `LFunction χ` from a `‖L χ c‖`-free bound on `B`

L-analogue of ζ §1444.  Using §1538 (`‖L χ c‖ ≥ 1/(3+1/(σ-1))`), `Mζ/‖L χ c‖ ≤ Mζ·(3+1/(σ-1))`, so §1509's
`hBexplicit` holds whenever `B ≥ 4·log(Mζ·(3+1/(σ-1)))/(R/2)` — an inequality with NO `‖L χ c‖`, ready for the
scalar constant choice (`Mζ ≍ |γ|`, `σ-1 = 1/(2K) ≍ 1/log|γ|`, RHS `≍ log|γ|`).  The `L χ c ≠ 0` positivity comes
from nonvanishing on `Re ≥ 1` (`LFunction_ne_zero_of_one_le_re`, χ≠1).  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem hBexplicit_of_ge_LFunction {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c : ℂ} {Mζ R B : ℝ} (hc : 1 < c.re) (hMζ : 0 < Mζ) (hR : 0 < R)
    (hB : 4 * Real.log (Mζ * (3 + 1 / (c.re - 1))) / (R / 2) ≤ B) :
    4 * Real.log (Mζ / ‖LFunction χ c‖) / (R / 2) ≤ B := by
  have hlb := norm_LFunction_ge_explicit χ hc
  have hζpos : 0 < ‖LFunction χ c‖ :=
    norm_pos_iff.mpr (LFunction_ne_zero_of_one_le_re χ (Or.inl hχ) hc.le)
  have hcre1 : 0 < c.re - 1 := by linarith
  have hXpos : 0 < 3 + 1 / (c.re - 1) := by positivity
  have hdiv : Mζ / ‖LFunction χ c‖ ≤ Mζ * (3 + 1 / (c.re - 1)) := by
    rw [div_le_iff₀ hζpos]
    have hmul : Mζ * (3 + 1 / (c.re - 1)) * (3 + 1 / (c.re - 1))⁻¹
        ≤ Mζ * (3 + 1 / (c.re - 1)) * ‖LFunction χ c‖ :=
      mul_le_mul_of_nonneg_left hlb (by positivity)
    rwa [mul_assoc, mul_inv_cancel₀ (ne_of_gt hXpos), mul_one] at hmul
  have hlog : Real.log (Mζ / ‖LFunction χ c‖) ≤ Real.log (Mζ * (3 + 1 / (c.re - 1))) :=
    Real.log_le_log (by positivity) hdiv
  have hRhalf : 0 < R / 2 := by linarith
  calc 4 * Real.log (Mζ / ‖LFunction χ c‖) / (R / 2)
      ≤ 4 * Real.log (Mζ * (3 + 1 / (c.re - 1))) / (R / 2) := by gcongr
    _ ≤ B := hB

#print axioms hBexplicit_of_ge_LFunction
