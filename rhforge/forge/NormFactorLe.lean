import Mathlib.Analysis.Normed.Field.Basic
import Mathlib.Analysis.SpecialFunctions.Complex.Circle

/-!
# §1373 — factor-norm bound: `‖g‖ ≤ Mζ / Plow` from `ζ = P·g`

The per-point arithmetic core of the factor sup-bound.  If `ζ = P·g` with `‖ζ‖ ≤ Mζ` and `‖P‖ ≥ Plow > 0`,
then the nonvanishing factor obeys `‖g‖ ≤ Mζ / Plow`.  Applied pointwise on the sphere (with `Mζ` from §1130 and
`Plow` a lower bound on the zero-product), then lifted to the disk by max-modulus §1370, this supplies the `M'`
that §1369 turns into §1355's `hgB`.  EXIT0 / std-axioms.
-/

theorem norm_factor_le {ζ P g : ℂ} {Mζ Plow : ℝ} (hfac : ζ = P * g) (hPpos : 0 < Plow)
    (hP : Plow ≤ ‖P‖) (hζ : ‖ζ‖ ≤ Mζ) : ‖g‖ ≤ Mζ / Plow := by
  have hnorm : ‖ζ‖ = ‖P‖ * ‖g‖ := by rw [hfac, norm_mul]
  rw [hnorm] at hζ
  rw [le_div_iff₀ hPpos]
  calc ‖g‖ * Plow ≤ ‖g‖ * ‖P‖ := mul_le_mul_of_nonneg_left hP (norm_nonneg g)
    _ = ‖P‖ * ‖g‖ := by ring
    _ ≤ Mζ := hζ

#print axioms norm_factor_le
