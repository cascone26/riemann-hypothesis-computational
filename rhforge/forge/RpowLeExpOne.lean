import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# §1292 — `N^a ≤ e` when `a·log N ≤ 1` (the `N^{1-σ}` collapse lemma)

The elementary inequality driving the collapse of §1291's growth bound: with `N ≍ |t|` and `σ ≥ 1 − 1/log N`,
the main terms `N^{1-σ}` are bounded by a constant, because `N^{1-σ} = exp((1-σ)·log N) ≤ exp 1 = e` once
`(1-σ)·log N ≤ 1`.  This is what turns the polynomial-looking `N^{1-σ}/(1-σ)` into log-type growth in `|t|`.
EXIT0 / std-axioms.
-/

theorem rpow_le_exp_one_of_mul_log_le {N : ℝ} (hN : 1 ≤ N) {a : ℝ}
    (h : a * Real.log N ≤ 1) : N ^ a ≤ Real.exp 1 := by
  rw [Real.rpow_def_of_pos (lt_of_lt_of_le one_pos hN)]
  exact Real.exp_le_exp.mpr (by rw [mul_comm]; exact h)

#print axioms rpow_le_exp_one_of_mul_log_le
