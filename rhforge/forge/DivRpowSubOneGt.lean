import ExpSubOneLeMulExp
import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# §1304 — `(N^{1-σ}-1)/(1-σ) ≤ log N` for `σ > 1` (σ>1 side of the two-sided collapse)

The mirror of §1294 for the `σ > 1` half of §1303's bound.  For `σ > 1`, `N^{1-σ} = e^{-x}` with
`x = (σ-1)·log N ≥ 0`, and `(N^{1-σ}-1)/(1-σ) = (1-e^{-x})/x · log N ≤ log N` since `1 - e^{-x} ≤ x`
(`Real.add_one_le_exp`).  Cleaner than the `σ<1` case — no `e` factor, no regime hypothesis.  Together with
§1295 this bounds §1303's first main term by `≲ log N` on BOTH sides of `σ = 1`, uniform over a BC disk.
EXIT0 / std-axioms.
-/

theorem div_rpow_sub_one_le_gt {N : ℝ} (hN : 1 ≤ N) {σ : ℝ} (hσ : 1 < σ) :
    (N ^ (1 - σ) - 1) / (1 - σ) ≤ Real.log N := by
  have hNpos : (0 : ℝ) < N := lt_of_lt_of_le one_pos hN
  have hL : 0 ≤ Real.log N := Real.log_nonneg hN
  set x := (σ - 1) * Real.log N with hx
  have hNpow : N ^ (1 - σ) = Real.exp (-x) := by
    rw [Real.rpow_def_of_pos hNpos, hx]; congr 1; ring
  have hexp : 1 - Real.exp (-x) ≤ x := by
    have h := Real.add_one_le_exp (-x); linarith
  have hid : Real.log N * (1 - σ) = -x := by rw [hx]; ring
  rw [hNpow, div_le_iff_of_neg (by linarith : (1 : ℝ) - σ < 0)]
  linarith [hexp, hid]

#print axioms div_rpow_sub_one_le_gt
