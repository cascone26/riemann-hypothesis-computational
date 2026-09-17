import DivRpowSubOneLe
import RpowLeExpOne

/-!
# §1295 — First main term collapses to `e·log N` in the collapse regime

Combining §1294 (`(N^{1-σ}-1)/(1-σ) ≤ N^{1-σ}·log N`) with §1292 (`N^{1-σ} ≤ e` when `(1-σ)·log N ≤ 1`):
in the regime `(1-σ)·log N ≤ 1` (i.e. `σ ≥ 1 − 1/log N`, the relevant approach to `σ = 1`), §1291's first main
term is bounded by `e·log N` — pure log-type growth.  With `N ≍ |t|` this is `≲ log|t|`.  EXIT0 / std-axioms.
-/

theorem main_term_le_exp_log {N : ℝ} (hN : 1 ≤ N) {σ : ℝ} (hσ : σ < 1)
    (hreg : (1 - σ) * Real.log N ≤ 1) :
    (N ^ (1 - σ) - 1) / (1 - σ) ≤ Real.exp 1 * Real.log N := by
  have hL : 0 ≤ Real.log N := Real.log_nonneg hN
  refine le_trans (div_rpow_sub_one_le hN hσ) ?_
  exact mul_le_mul_of_nonneg_right (rpow_le_exp_one_of_mul_log_le hN hreg) hL

#print axioms main_term_le_exp_log
