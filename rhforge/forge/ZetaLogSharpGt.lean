import ZetaGrowthTwoSided
import DivRpowSubOneGt
import Term4SharpBound

/-!
# §1305 — Sharp ζ log bound for `σ > 1`: `‖ζ(s)‖ ≤ log N + 2 + 1/‖s-1‖ + 1/σ`

The `σ > 1` companion of §1298, using the two-sided bound §1303 and the `σ>1` collapse §1304.  For `σ > 1`,
`N^{1-σ} ≤ 1` (since `N ≥ 1`, `1-σ ≤ 0`), so every remaining term drops with no regime hypothesis:
term 1 `≤ log N` (§1304), term 3 `N^{1-σ}/‖s-1‖ ≤ 1/‖s-1‖`, term 4 `≤ 1 + N^{1-σ}/σ ≤ 1 + 1/σ` (§1297).
Requires `|Im s| ≤ N`.  Together with §1298 this covers the whole `σ ≠ 1` range for a BC disk.  EXIT0 / std-axioms.
-/

theorem norm_riemannZeta_le_log_gt {N : ℕ} (hN : 1 ≤ N) {s : ℂ}
    (hs1 : 1 < s.re) (hsim : |s.im| ≤ (N : ℝ)) :
    ‖riemannZeta s‖ ≤ Real.log N + 2 + 1 / ‖s - 1‖ + 1 / s.re := by
  have hs0 : 0 < s.re := lt_trans one_pos hs1
  have hs1ne : s.re ≠ 1 := ne_of_gt hs1
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  have hpow_le1 : (N : ℝ) ^ (1 - s.re) ≤ 1 :=
    Real.rpow_le_one_of_one_le_of_nonpos hNR (by linarith)
  have hd : (0 : ℝ) ≤ ‖s - 1‖ := norm_nonneg _
  have hT1 : ((N : ℝ) ^ (1 - s.re) - 1) / (1 - s.re) ≤ Real.log N :=
    div_rpow_sub_one_le_gt hNR hs1
  have hT2 : (N : ℝ) ^ (1 - s.re) / ‖s - 1‖ ≤ 1 / ‖s - 1‖ :=
    div_le_div_of_nonneg_right hpow_le1 hd
  have hT3a : ‖s‖ * (N : ℝ) ^ (-s.re) / s.re ≤ 1 + (N : ℝ) ^ (1 - s.re) / s.re :=
    norm_mul_rpow_neg_div_le hs0 hNR hsim
  have hT3b : (N : ℝ) ^ (1 - s.re) / s.re ≤ 1 / s.re :=
    div_le_div_of_nonneg_right hpow_le1 hs0.le
  have hstrip : ‖riemannZeta s‖ ≤ ((N : ℝ) ^ (1 - s.re) - 1) / (1 - s.re) + 1
      + (N : ℝ) ^ (1 - s.re) / ‖s - 1‖ + ‖s‖ * ((N : ℝ) ^ (-s.re) / s.re) :=
    norm_riemannZeta_le_offline hN hs0 hs1ne
  rw [← mul_div_assoc] at hstrip
  linarith [hstrip, hT1, hT2, hT3a, hT3b]

#print axioms norm_riemannZeta_le_log_gt
