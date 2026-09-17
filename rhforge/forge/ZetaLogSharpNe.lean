import ZetaLogSharp
import ZetaLogSharpGt

/-!
# §1306 — Unified sharp ζ log bound for `σ ≠ 1` (both sides): `‖ζ(s)‖ ≤ e·log N + 2 + e/‖s-1‖ + e/σ`

Merging §1298 (`σ < 1`) and §1305 (`σ > 1`) into one statement covering the whole `σ ≠ 1` range.  On the `σ > 1`
side §1305's bound `log N + 2 + 1/‖s-1‖ + 1/σ` is dominated by the `σ < 1` shape since `e ≥ 1`
(`log N ≤ e·log N`, `1/x ≤ e/x`).  Result, for `0 < σ`, `σ ≠ 1`, `|Im s| ≤ N`, region `(1-σ)·log N ≤ 1`:

  `‖ζ(s)‖ ≤ e·log N + 2 + e/‖s-1‖ + e/σ`.

A single bound valid across a Borel–Carathéodory disk that straddles `σ = 1` — the uniform-disk ingredient
feeding `deriv_bound_of_re_le_general`.  EXIT0 / std-axioms.
-/

theorem norm_riemannZeta_le_log_ne {N : ℕ} (hN : 1 ≤ N) {s : ℂ}
    (hs0 : 0 < s.re) (hs1 : s.re ≠ 1) (hsim : |s.im| ≤ (N : ℝ))
    (hreg : (1 - s.re) * Real.log N ≤ 1) :
    ‖riemannZeta s‖ ≤ Real.exp 1 * Real.log N + 2 + Real.exp 1 / ‖s - 1‖ + Real.exp 1 / s.re := by
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  have he1 : (1 : ℝ) ≤ Real.exp 1 := Real.one_le_exp_iff.mpr zero_le_one
  rcases lt_or_gt_of_ne hs1 with hlt | hgt
  · exact norm_riemannZeta_le_log_sharp hN hs0 hlt hreg hsim
  · have h := norm_riemannZeta_le_log_gt hN hgt hsim
    have hsne1 : s ≠ 1 := fun he => by rw [he] at hgt; simp at hgt
    have hdpos : (0 : ℝ) < ‖s - 1‖ := norm_pos_iff.mpr (sub_ne_zero.mpr hsne1)
    have hlogN : (0 : ℝ) ≤ Real.log N := Real.log_nonneg hNR
    have b1 : Real.log N ≤ Real.exp 1 * Real.log N := le_mul_of_one_le_left hlogN he1
    have b2 : 1 / ‖s - 1‖ ≤ Real.exp 1 / ‖s - 1‖ := div_le_div_of_nonneg_right he1 hdpos.le
    have b3 : 1 / s.re ≤ Real.exp 1 / s.re := div_le_div_of_nonneg_right he1 hs0.le
    linarith [h, b1, b2, b3]

#print axioms norm_riemannZeta_le_log_ne
