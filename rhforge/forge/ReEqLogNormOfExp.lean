import Mathlib.Analysis.SpecialFunctions.Complex.Log

/-!
# §1320 — Any log-branch has real part `log‖h‖`: `exp w = v ⟹ Re w = log‖v‖`

The real part of a holomorphic branch of `log h` is branch-INDEPENDENT: if `exp w = v` then
`‖v‖ = ‖exp w‖ = exp(Re w)`, so `Re w = log‖v‖`.  This is why the Borel–Carathéodory bound (which only sees
`Re f = log‖h‖ ≤ M`) is well-posed even where the principal `Complex.log` branch fails (slitPlane) — ANY branch
gives the same `Re`.  It isolates the sole remaining ζ-analytic gap to pure branch EXISTENCE (a holomorphic `f`
with `exp∘f = ζ/P` on the disk — the simply-connected log), the value/derivative facts being branch-free.
EXIT0 / std-axioms.
-/

open Complex

theorem re_eq_log_norm_of_exp {w v : ℂ} (h : Complex.exp w = v) : w.re = Real.log ‖v‖ := by
  rw [← h, Complex.norm_exp, Real.log_exp]

#print axioms re_eq_log_norm_of_exp
