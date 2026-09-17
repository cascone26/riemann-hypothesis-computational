import Mathlib.Analysis.SpecialFunctions.Complex.Log

/-!
# §1319 — `Re(log w) ≤ log M'` from `‖w‖ ≤ M'` (the `Re(log h) ≤ M` input to §1318)

`Re(Complex.log w) = log‖w‖` (`Complex.log_re`), so a norm bound `‖w‖ ≤ M'` gives `Re(log w) ≤ log M'` by
monotonicity of `Real.log`.  Applied pointwise to `w = (ζ/P)(z)` on the disk (bounded by the ball bound §1130),
this is exactly the `Re(f z) ≤ M` hypothesis (`M = log M'`) that Borel–Carathéodory §1318 consumes.
EXIT0 / std-axioms.
-/

open Complex

theorem re_log_le_of_norm_le {w : ℂ} {M' : ℝ} (hw : w ≠ 0) (hb : ‖w‖ ≤ M') :
    (Complex.log w).re ≤ Real.log M' := by
  rw [Complex.log_re]
  exact Real.log_le_log (norm_pos_iff.mpr hw) hb

#print axioms re_log_le_of_norm_le
