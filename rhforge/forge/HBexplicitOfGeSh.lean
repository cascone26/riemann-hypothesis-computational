import NormZetaGeExplicit
import Mathlib.NumberTheory.LSeries.Dirichlet

/-! [S3020 SHARP COPY: Borel–Carathéodory constant 4 → 2 via SharpBC (S3018); names suffixed _sh]

# §1444 — discharge the BC scalar `hBexplicit` from a `‖ζ c‖`-free bound on `B`

Using §1443 (`‖ζ c‖ ≥ 1/(3+1/(σ-1))`), `Mζ/‖ζ c‖ ≤ Mζ·(3+1/(σ-1))`, so
`4·log(Mζ/‖ζ c‖)/(R/2) ≤ 4·log(Mζ·(3+1/(σ-1)))/(R/2)`.  Hence §1437's `hBexplicit` holds whenever
`B ≥ 4·log(Mζ·(3+1/(σ-1)))/(R/2)` — an inequality with NO `‖ζ c‖`, ready for the scalar constant choice
(`Mζ ≍ |γ|`, `σ-1 = 1/(2K) ≍ 1/log|γ|`, so the RHS `≍ log|γ|`).  EXIT0 / std-axioms.
-/

open Complex

theorem hBexplicit_of_ge_sh {c : ℂ} {Mζ R B : ℝ} (hc : 1 < c.re) (hMζ : 0 < Mζ) (hR : 0 < R)
    (hB : 2 * Real.log (Mζ * (3 + 1 / (c.re - 1))) / (R / 2) ≤ B) :
    2 * Real.log (Mζ / ‖riemannZeta c‖) / (R / 2) ≤ B := by
  have hlb := norm_riemannZeta_ge_explicit hc
  have hζpos : 0 < ‖riemannZeta c‖ := norm_pos_iff.mpr (riemannZeta_ne_zero_of_one_lt_re hc)
  have hcre1 : 0 < c.re - 1 := by linarith
  have hXpos : 0 < 3 + 1 / (c.re - 1) := by positivity
  have hdiv : Mζ / ‖riemannZeta c‖ ≤ Mζ * (3 + 1 / (c.re - 1)) := by
    rw [div_le_iff₀ hζpos]
    have hmul : Mζ * (3 + 1 / (c.re - 1)) * (3 + 1 / (c.re - 1))⁻¹
        ≤ Mζ * (3 + 1 / (c.re - 1)) * ‖riemannZeta c‖ :=
      mul_le_mul_of_nonneg_left hlb (by positivity)
    rwa [mul_assoc, mul_inv_cancel₀ (ne_of_gt hXpos), mul_one] at hmul
  have hlog : Real.log (Mζ / ‖riemannZeta c‖) ≤ Real.log (Mζ * (3 + 1 / (c.re - 1))) :=
    Real.log_le_log (by positivity) hdiv
  have hRhalf : 0 < R / 2 := by linarith
  calc 2 * Real.log (Mζ / ‖riemannZeta c‖) / (R / 2)
      ≤ 2 * Real.log (Mζ * (3 + 1 / (c.re - 1))) / (R / 2) := by gcongr
    _ ≤ B := hB

#print axioms hBexplicit_of_ge_sh
