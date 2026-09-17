import ZetaStripGrowth
import MainTermLogBound
import Term4SharpBound

/-!
# §1298 — Sharp ζ log bound (corrected tail): `‖ζ(s)‖ ≤ e·log N + 2 + e/‖s-1‖ + e/σ`

The honest log-type bound, using the paired tail §1297 instead of §1296's lossy `N^{-σ}≤1`.  For `0<σ<1`,
`(1-σ)·log N ≤ 1`, and `|Im s| ≤ N`:
* term 1 `≤ e·log N` (§1295),
* term 3 `N^{1-σ}/‖s-1‖ ≤ e/‖s-1‖` (§1292),
* term 4 `‖s‖·N^{-σ}/σ ≤ 1 + N^{1-σ}/σ ≤ 1 + e/σ` (§1297 + §1292),

giving `‖ζ(s)‖ ≤ e·log N + 2 + e/‖s-1‖ + e/σ`.  At `N = ⌈|Im s|⌉` (so `|Im s| ≤ N`), with `‖s-1‖ ≥ |Im s|`,
this is genuinely `‖ζ(σ+it)‖ ≲ log|t|` — no spurious `|t|`-growth.  EXIT0 / std-axioms.
-/

theorem norm_riemannZeta_le_log_sharp {N : ℕ} (hN : 1 ≤ N) {s : ℂ}
    (hs0 : 0 < s.re) (hs1 : s.re < 1) (hreg : (1 - s.re) * Real.log N ≤ 1)
    (hsim : |s.im| ≤ (N : ℝ)) :
    ‖riemannZeta s‖ ≤ Real.exp 1 * Real.log N + 2 + Real.exp 1 / ‖s - 1‖ + Real.exp 1 / s.re := by
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  have he : (N : ℝ) ^ (1 - s.re) ≤ Real.exp 1 := rpow_le_exp_one_of_mul_log_le hNR hreg
  have hd : (0 : ℝ) ≤ ‖s - 1‖ := norm_nonneg _
  have hT1 : ((N : ℝ) ^ (1 - s.re) - 1) / (1 - s.re) ≤ Real.exp 1 * Real.log N :=
    main_term_le_exp_log hNR hs1 hreg
  have hT2 : (N : ℝ) ^ (1 - s.re) / ‖s - 1‖ ≤ Real.exp 1 / ‖s - 1‖ :=
    div_le_div_of_nonneg_right he hd
  have hT3a : ‖s‖ * (N : ℝ) ^ (-s.re) / s.re ≤ 1 + (N : ℝ) ^ (1 - s.re) / s.re :=
    norm_mul_rpow_neg_div_le hs0 hNR hsim
  have hT3b : (N : ℝ) ^ (1 - s.re) / s.re ≤ Real.exp 1 / s.re :=
    div_le_div_of_nonneg_right he hs0.le
  have hstrip : ‖riemannZeta s‖ ≤ ((N : ℝ) ^ (1 - s.re) - 1) / (1 - s.re) + 1
      + (N : ℝ) ^ (1 - s.re) / ‖s - 1‖ + ‖s‖ * ((N : ℝ) ^ (-s.re) / s.re) :=
    norm_riemannZeta_le_strip_EM hN hs0 hs1
  rw [← mul_div_assoc] at hstrip
  linarith [hstrip, hT1, hT2, hT3a, hT3b]

#print axioms norm_riemannZeta_le_log_sharp
