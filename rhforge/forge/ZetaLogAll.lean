import ZetaLogSharp
import ZetaLogEqOne
import ZetaLogSharpGt

/-!
# §1309 ★ Fully uniform ζ log bound on a BC disk: `‖ζ(s)‖ ≤ e·log N + 3 + e/‖s-1‖ + e/σ`

Merging all three σ-regimes — §1298 (`σ<1`), §1308 (`σ=1`), §1305 (`σ>1`) — into a single bound valid for
every `0 < σ = Re s`, `s ≠ 1`, `|Im s| ≤ N`, with the region hypothesis `(1-σ)·log N ≤ 1` (which is free for
`σ ≥ 1`).  Since `e ≥ 1` dominates the `σ ≥ 1` constants, the common form is

  `‖ζ(s)‖ ≤ e·log N + 3 + e/‖s-1‖ + e/σ`.

This is THE uniform bound over a Borel–Carathéodory disk centred near `1 + it₀` (straddling `σ = 1`), the exact
`Re(Log ζ) ≤ M` input for `deriv_bound_of_re_le_general`.  EXIT0 / std-axioms.
-/

theorem norm_riemannZeta_le_log_all {N : ℕ} (hN : 1 ≤ N) {s : ℂ}
    (hs0 : 0 < s.re) (hsne : s ≠ 1) (hsim : |s.im| ≤ (N : ℝ))
    (hreg : (1 - s.re) * Real.log N ≤ 1) :
    ‖riemannZeta s‖ ≤ Real.exp 1 * Real.log N + 3 + Real.exp 1 / ‖s - 1‖ + Real.exp 1 / s.re := by
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  have he1 : (1 : ℝ) ≤ Real.exp 1 := Real.one_le_exp_iff.mpr zero_le_one
  have hlogN : (0 : ℝ) ≤ Real.log N := Real.log_nonneg hNR
  have hdpos : (0 : ℝ) < ‖s - 1‖ := norm_pos_iff.mpr (sub_ne_zero.mpr hsne)
  have hediv : (0 : ℝ) ≤ Real.exp 1 / s.re := (div_pos (Real.exp_pos 1) hs0).le
  have hb1 : Real.log N ≤ Real.exp 1 * Real.log N := le_mul_of_one_le_left hlogN he1
  have hb2 : 1 / ‖s - 1‖ ≤ Real.exp 1 / ‖s - 1‖ := div_le_div_of_nonneg_right he1 hdpos.le
  have hb3 : 1 / s.re ≤ Real.exp 1 / s.re := div_le_div_of_nonneg_right he1 hs0.le
  rcases lt_trichotomy s.re 1 with hlt | heq | hgt
  · have h := norm_riemannZeta_le_log_sharp hN hs0 hlt hreg hsim
    linarith [h]
  · have h := norm_riemannZeta_le_log_eq_one hN heq hsne hsim
    linarith [h, hb1, hb2, hediv]
  · have h := norm_riemannZeta_le_log_gt hN hgt hsim
    linarith [h, hb1, hb2, hb3]

#print axioms norm_riemannZeta_le_log_all
