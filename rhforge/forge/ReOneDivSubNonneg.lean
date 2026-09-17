import Mathlib.Analysis.SpecialFunctions.Complex.Circle

/-!
# §1396 — `0 ≤ Re(1/(z-ρ))` when `Re ρ ≤ Re z`

`Re(1/(z-ρ)) = (Re z - Re ρ)/|z-ρ|² ≥ 0` whenever `Re ρ ≤ Re z`.  This discharges §1395's `hre` at the growth point
`z = σ+2iγ` (`σ > 1 ≥ Re ρ` for every ζ-zero `ρ`), feeding the drop-zeros bound `Re(-logDeriv ζ) ≤ -Re(logDeriv
g_full)` for `hgrowth`.  EXIT0 / std-axioms.
-/

open Complex

theorem re_one_div_sub_nonneg {z ρ : ℂ} (h : ρ.re ≤ z.re) : 0 ≤ (1 / (z - ρ)).re := by
  rw [one_div, Complex.inv_re, Complex.sub_re]
  apply div_nonneg
  · linarith
  · exact Complex.normSq_nonneg _

#print axioms re_one_div_sub_nonneg
