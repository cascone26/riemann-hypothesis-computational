import ExpSubOneLeMulExp
import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# §1294 — `(N^{1-σ}-1)/(1-σ) ≤ N^{1-σ}·log N` (collapsing §1291's first main term)

The divided-difference form `(e^y-1)/y ≤ e^y` (§1293) at `y = (1-σ)·log N` turns §1291's main-sum bound term
`(N^{1-σ}-1)/(1-σ)` into `N^{1-σ}·log N` — genuine log-type growth once `N^{1-σ}` is bounded (§1292).  Valid for
`N ≥ 1`, `σ < 1`.  EXIT0 / std-axioms.
-/

theorem div_rpow_sub_one_le {N : ℝ} (hN : 1 ≤ N) {σ : ℝ} (hσ : σ < 1) :
    (N ^ (1 - σ) - 1) / (1 - σ) ≤ N ^ (1 - σ) * Real.log N := by
  set a := 1 - σ with ha
  have ha0 : 0 < a := by rw [ha]; linarith
  have hL : 0 ≤ Real.log N := Real.log_nonneg hN
  have hNa : N ^ a = Real.exp (Real.log N * a) :=
    Real.rpow_def_of_pos (lt_of_lt_of_le one_pos hN) a
  have hy : 0 ≤ Real.log N * a := mul_nonneg hL (le_of_lt ha0)
  have hexp := exp_sub_one_le_mul_exp hy
  rw [← hNa] at hexp
  rw [div_le_iff₀ ha0]
  nlinarith [hexp]

#print axioms div_rpow_sub_one_le
