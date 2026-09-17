import Mathlib.Analysis.SpecialFunctions.Log.Basic

/-!
# §1456 — bound §1454's `hBexp` RHS by `32·log((3|γ|+10)(3+2K))`

With `σK = 1+1/(2K) ∈ (1,2]` (`K ≥ 1/2`) and `|γ| > 1`: `1/(σK-1) = 2K`, the inner factor
`3 + 1/(|γ|-1/4) + (σK+|γ|+1/4)/(σK-1/4) ≤ 3|γ|+10`, and `4·(…)/(1/4/2) = 32·(…)`.  So the giant `hBexp` RHS is
`≤ 32·log((3|γ|+10)(3+2K))` — explicit and linear in `|γ|`, ready for the §1445 threshold.  EXIT0 / std-axioms.
-/

theorem hBexp_rhs_le {K γ : ℝ} (hK : (1 : ℝ) / 2 ≤ K) (hγ : 1 < |γ|) :
    4 * Real.log ((3 + 1 / (|γ| - 1 / 4)
        + ((1 + 1 / (2 * K)) + |γ| + 1 / 4) / ((1 + 1 / (2 * K)) - 1 / 4))
        * (3 + 1 / ((1 + 1 / (2 * K)) - 1))) / (1 / 4 / 2)
      ≤ 32 * Real.log ((3 * |γ| + 10) * (3 + 2 * K)) := by
  have hK0 : 0 < K := by linarith
  have hσ1 : 1 < 1 + 1 / (2 * K) := by
    have h : 0 < 1 / (2 * K) := by positivity
    linarith
  have hσ2 : 1 + 1 / (2 * K) ≤ 2 := by
    have : 1 / (2 * K) ≤ 1 := by rw [div_le_one (by positivity)]; linarith
    linarith
  have hden_pos : 0 < (1 + 1 / (2 * K)) - 1 / 4 := by linarith
  have hγpos : 0 < |γ| - 1 / 4 := by linarith
  have hinv : 1 / ((1 + 1 / (2 * K)) - 1) = 2 * K := by
    rw [show (1 + 1 / (2 * K)) - 1 = 1 / (2 * K) by ring, one_div_one_div]
  have ha : 1 / (|γ| - 1 / 4) ≤ 4 / 3 := by rw [div_le_iff₀ hγpos]; nlinarith [hγ]
  have hb : ((1 + 1 / (2 * K)) + |γ| + 1 / 4) / ((1 + 1 / (2 * K)) - 1 / 4) ≤ (4 / 3) * (|γ| + 9 / 4) := by
    rw [div_le_iff₀ hden_pos]
    nlinarith [hσ2, hden_pos, abs_nonneg γ, hγ,
      mul_nonneg (show (0 : ℝ) ≤ (4 / 3) * (|γ| + 9 / 4) by positivity)
        (show (0 : ℝ) ≤ ((1 + 1 / (2 * K)) - 1 / 4) - 3 / 4 by linarith)]
  have hinner : 3 + 1 / (|γ| - 1 / 4)
      + ((1 + 1 / (2 * K)) + |γ| + 1 / 4) / ((1 + 1 / (2 * K)) - 1 / 4) ≤ 3 * |γ| + 10 := by
    nlinarith [ha, hb, abs_nonneg γ, hγ]
  have hInnerPos : 0 < 3 + 1 / (|γ| - 1 / 4)
      + ((1 + 1 / (2 * K)) + |γ| + 1 / 4) / ((1 + 1 / (2 * K)) - 1 / 4) := by positivity
  have hKfac : 0 < 3 + 2 * K := by linarith
  have hlog : Real.log ((3 + 1 / (|γ| - 1 / 4)
        + ((1 + 1 / (2 * K)) + |γ| + 1 / 4) / ((1 + 1 / (2 * K)) - 1 / 4)) * (3 + 2 * K))
      ≤ Real.log ((3 * |γ| + 10) * (3 + 2 * K)) :=
    Real.log_le_log (by positivity) (by nlinarith [hinner, hKfac, hInnerPos])
  rw [hinv]
  have h32 : 4 * Real.log ((3 + 1 / (|γ| - 1 / 4)
        + ((1 + 1 / (2 * K)) + |γ| + 1 / 4) / ((1 + 1 / (2 * K)) - 1 / 4)) * (3 + 2 * K)) / (1 / 4 / 2)
      = 32 * Real.log ((3 + 1 / (|γ| - 1 / 4)
        + ((1 + 1 / (2 * K)) + |γ| + 1 / 4) / ((1 + 1 / (2 * K)) - 1 / 4)) * (3 + 2 * K)) := by
    ring
  rw [h32]
  linarith [hlog]

#print axioms hBexp_rhs_le
