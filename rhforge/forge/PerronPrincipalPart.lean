import PerronPoleRemovable
import Mathlib.Analysis.Calculus.DSlope

set_option maxHeartbeats 1000000

/-!
# §1614 — the principal-part split: `dslope Ĥ 1 s = F(s) - (x²/2)/(s-1)` (`s ≠ 1`)

With `Ĥ = update ((s-1)·F) 1 (x²/2)` (analytic at 1, §1613), its slope at `1` recovers the holomorphic remainder of
the Laurent split: `dslope Ĥ 1 s = (Ĥ s - Ĥ 1)/(s-1) = ((s-1)F(s) - x²/2)/(s-1) = F(s) - (x²/2)/(s-1)`.  Equivalently
`F(s) = (x²/2)/(s-1) + dslope Ĥ 1 s`: the Perron integrand equals its explicit `s=1` pole term plus a piece that is
analytic at `1` — the decomposition the contour shift needs to deposit exactly the residue `x²/2`.  EXIT0 / std.
-/

open Complex

theorem perron_dslope_eq {x : ℝ} {s : ℂ} (hs1 : s ≠ 1) :
    dslope
        (Function.update
          (fun z => (z - 1) * (-logDeriv riemannZeta z * (x : ℂ) ^ (z + 1) / (z * (z + 1))))
          1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1 s
      = -logDeriv riemannZeta s * (x : ℂ) ^ (s + 1) / (s * (s + 1)) - (x : ℂ) ^ (2 : ℂ) / 2 / (s - 1) := by
  have hsub : s - 1 ≠ 0 := sub_ne_zero.mpr hs1
  rw [dslope_of_ne _ hs1, slope_def_field, Function.update_of_ne hs1, Function.update_self]
  generalize -logDeriv riemannZeta s * (x : ℂ) ^ (s + 1) / (s * (s + 1)) = A
  field_simp

#print axioms perron_dslope_eq
