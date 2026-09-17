import Mathlib.Analysis.Complex.Norm

/-!
# §1317 — `‖g‖ ≤ B ⟹ -B ≤ Re g` (the final connector: BC bound → §1314 input)

The Borel–Carathéodory derivative bound gives `‖logDeriv h‖ ≤ B`; §1314 needs `-B ≤ Re(logDeriv h)`.  This is the
one-line bridge: `Re g ≥ -|Re g| ≥ -‖g‖ ≥ -B`.  With it, the sole remaining analytic content of the whole
quantitative zero-free region is producing the norm bound `‖logDeriv(ζ/P)‖ ≤ B` itself (deriv_bound on `log h`).
EXIT0 / std-axioms.
-/

open Complex

theorem neg_le_re_of_norm_le {g : ℂ} {B : ℝ} (h : ‖g‖ ≤ B) : -B ≤ g.re := by
  have h1 : -‖g‖ ≤ g.re := (abs_le.mp (Complex.abs_re_le_norm g)).1
  linarith

#print axioms neg_le_re_of_norm_le
