import TwoModeParseval

/-!
# §1736 — The archimedean metric is positive definite (two-mode)

Immediate from §1735's Parseval value `(c²+d²)(L/2)`: on a window of positive length the only two-mode
combination with zero archimedean self-energy is the trivial one.

`0 < L` and `∫₀ᴸ (c·sin(a u) + d·sin(b u))² du = 0  ⟹  c = 0 ∧ d = 0`.

This is the genuine positive-definiteness (not just nonnegativity) of the `L²` metric on the interval
eigenbasis — the property the Rosati/critical-circle argument (§1729) needs of its inner product.  EXIT0.
-/

open Real intervalIntegral MeasureTheory

/-- **Two-mode positive-definiteness.**  Zero self-energy forces zero coefficients. -/
theorem two_mode_metric_pos_def (a b L c d : ℝ) (k l : ℤ)
    (hL : 0 < L) (ha : a ≠ 0) (hb : b ≠ 0) (hd : a - b ≠ 0) (hs : a + b ≠ 0)
    (hLa : a * L = (k : ℝ) * π) (hLb : b * L = (l : ℝ) * π)
    (h0 : (∫ u in (0:ℝ)..L, (c * Real.sin (a * u) + d * Real.sin (b * u)) ^ 2) = 0) :
    c = 0 ∧ d = 0 := by
  rw [two_mode_parseval a b L c d k l ha hb hd hs hLa hLb] at h0
  have hL2 : 0 < L / 2 := by linarith
  have hcd : c ^ 2 + d ^ 2 = 0 := by
    rcases mul_eq_zero.mp h0 with h | h
    · exact h
    · exact absurd h (ne_of_gt hL2)
  exact ⟨by nlinarith [sq_nonneg c, sq_nonneg d], by nlinarith [sq_nonneg c, sq_nonneg d]⟩
