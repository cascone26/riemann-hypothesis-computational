import ZetaStripGrowth
import MainTermLogBound

/-!
# §1296 — ζ strip bound consolidated to log-type: no `N^{±σ}` powers, only `log N`

In the collapse regime `(1-σ)·log N ≤ 1` (i.e. `σ ≥ 1 − 1/log N`, `σ = Re s`), every `N`-power in §1291's bound
is discharged: the first main term `≤ e·log N` (§1295), `N^{1-σ} ≤ e` (§1292), `N^{-σ} ≤ 1`.  Result:

  `‖ζ(s)‖ ≤ e·log N + 1 + e/‖s-1‖ + ‖s‖/σ`   for `0 < σ < 1`, `(1-σ)·log N ≤ 1`.

Only `log N` survives as a growing quantity; the rest is controlled by `‖s‖, ‖s-1‖, σ`.  With `N ≍ |t|`
(so `log N ≍ log|t|`, `‖s-1‖ ≥ |t|`, `‖s‖ ≍ |t|`, `σ` near 1) this is `‖ζ(σ+it)‖ ≲ log|t|` — the sharp-region
input to §1269.  EXIT0 / std-axioms.
-/

theorem norm_riemannZeta_le_log {N : ℕ} (hN : 1 ≤ N) {s : ℂ}
    (hs0 : 0 < s.re) (hs1 : s.re < 1) (hreg : (1 - s.re) * Real.log N ≤ 1) :
    ‖riemannZeta s‖ ≤ Real.exp 1 * Real.log N + 1 + Real.exp 1 / ‖s - 1‖ + ‖s‖ / s.re := by
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  have hsne1 : s ≠ 1 := by intro h; rw [h] at hs1; simp at hs1
  have hd : (0 : ℝ) < ‖s - 1‖ := norm_pos_iff.mpr (sub_ne_zero.mpr hsne1)
  have hsn : (0 : ℝ) ≤ ‖s‖ := norm_nonneg _
  have he : (N : ℝ) ^ (1 - s.re) ≤ Real.exp 1 := rpow_le_exp_one_of_mul_log_le hNR hreg
  have hone : (N : ℝ) ^ (-s.re) ≤ 1 :=
    Real.rpow_le_one_of_one_le_of_nonpos hNR (neg_nonpos.mpr (le_of_lt hs0))
  have hT1 : ((N : ℝ) ^ (1 - s.re) - 1) / (1 - s.re) ≤ Real.exp 1 * Real.log N :=
    main_term_le_exp_log hNR hs1 hreg
  have hT2 : (N : ℝ) ^ (1 - s.re) / ‖s - 1‖ ≤ Real.exp 1 / ‖s - 1‖ :=
    div_le_div_of_nonneg_right he hd.le
  have hT3 : ‖s‖ * ((N : ℝ) ^ (-s.re) / s.re) ≤ ‖s‖ / s.re := by
    rw [← mul_div_assoc]
    refine div_le_div_of_nonneg_right ?_ hs0.le
    calc ‖s‖ * (N : ℝ) ^ (-s.re) ≤ ‖s‖ * 1 := mul_le_mul_of_nonneg_left hone hsn
      _ = ‖s‖ := mul_one _
  calc ‖riemannZeta s‖
      ≤ ((N : ℝ) ^ (1 - s.re) - 1) / (1 - s.re) + 1
          + (N : ℝ) ^ (1 - s.re) / ‖s - 1‖ + ‖s‖ * ((N : ℝ) ^ (-s.re) / s.re) :=
        norm_riemannZeta_le_strip_EM hN hs0 hs1
    _ ≤ Real.exp 1 * Real.log N + 1 + Real.exp 1 / ‖s - 1‖ + ‖s‖ / s.re := by
        gcongr

#print axioms norm_riemannZeta_le_log
