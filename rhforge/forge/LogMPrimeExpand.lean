import Mathlib.Analysis.SpecialFunctions.Log.Basic

/-!
# §1379 — `log M'` expansion: `log(Mζ/(R/2)^n) = log Mζ + n·log(2/R)`

The arithmetic making §1378's Borel–Carathéodory constant explicit and `≍ log|t|`.  With `M' = Mζ/(R/2)^n`,

  `log M' = log Mζ - n·log(R/2) = log Mζ + n·log(2/R)`.

Since `Mζ ≍ log|t|` (§1130/§1299 give `log Mζ ≍ log|t|`) and the half-disk zero count `n ≍ log|t|` (Jensen §1363),
`log M' ≍ log|t|`, so §1378's `B = 4(log M' - log‖g c‖)/(R/2)` is the `≍ log|t|` bound §1355 needs.  EXIT0 / std-axioms.
-/

theorem log_div_pow_half {Mζ R : ℝ} {n : ℕ} (hMζ : 0 < Mζ) (hR : 0 < R) :
    Real.log (Mζ / (R / 2) ^ n) = Real.log Mζ + n * Real.log (2 / R) := by
  have hR2 : (0 : ℝ) < R / 2 := by linarith
  rw [Real.log_div (ne_of_gt hMζ) (ne_of_gt (pow_pos hR2 n)), Real.log_pow,
    show (2 : ℝ) / R = (R / 2)⁻¹ from (inv_div R 2).symm, Real.log_inv]
  ring

#print axioms log_div_pow_half
