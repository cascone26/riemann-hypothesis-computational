import Mathlib.Analysis.SpecialFunctions.Log.Basic

set_option maxHeartbeats 1000000

/-!
# §1546 — bound §1545's L `hBexp` RHS by `32·log(N·(2|γ|+6)·(3+2K))`

L-analogue of ζ §1456, for the L inner factor `M'_L = 2 + N + (σK+|γ|+1/4)·N/(σK-1/4)` (modulus `N`, NO pole term).
With `σK = 1+1/(2K) ∈ (1,2]` (`K ≥ 1/2`) and `|γ| > 1`, `N ≥ 1`: `(σK+|γ|+1/4)/(σK-1/4) ≤ (4/3)|γ|+3`, so
`M'_L ≤ N·(2|γ|+6)`; and `4·(…)/(1/4/2) = 32·(…)`.  Hence the giant L `hBexp` RHS is `≤ 32·log(N·(2|γ|+6)·(3+2K))`
— explicit, linear in `|γ|` (with the fixed modulus `N`), ready for the §1445 threshold.  EXIT0 / std-axioms.
-/

theorem hBexp_rhs_le_LFunction {K γ : ℝ} {N : ℕ} (hK : (1 : ℝ) / 2 ≤ K) (hγ : 1 < |γ|)
    (hN1 : (1 : ℝ) ≤ (N : ℝ)) :
    4 * Real.log ((2 + (N : ℝ)
        + ((1 + 1 / (2 * K)) + |γ| + 1 / 4) * (N : ℝ) / ((1 + 1 / (2 * K)) - 1 / 4))
        * (3 + 1 / ((1 + 1 / (2 * K)) - 1))) / (1 / 4 / 2)
      ≤ 32 * Real.log ((N : ℝ) * (2 * |γ| + 6) * (3 + 2 * K)) := by
  have hK0 : 0 < K := by linarith
  have hN0 : (0 : ℝ) ≤ (N : ℝ) := by linarith
  have hσ1 : 1 < 1 + 1 / (2 * K) := by
    have h : 0 < 1 / (2 * K) := by positivity
    linarith
  have hσ2 : 1 + 1 / (2 * K) ≤ 2 := by
    have : 1 / (2 * K) ≤ 1 := by rw [div_le_one (by positivity)]; linarith
    linarith
  have hden_pos : 0 < (1 + 1 / (2 * K)) - 1 / 4 := by linarith
  have hγpos : 0 < |γ| := by linarith
  have hinv : 1 / ((1 + 1 / (2 * K)) - 1) = 2 * K := by
    rw [show (1 + 1 / (2 * K)) - 1 = 1 / (2 * K) by ring, one_div_one_div]
  -- ratio bound
  have hratio : ((1 + 1 / (2 * K)) + |γ| + 1 / 4) / ((1 + 1 / (2 * K)) - 1 / 4)
      ≤ (4 / 3) * |γ| + 3 := by
    rw [div_le_iff₀ hden_pos]
    nlinarith [hσ1, hσ2, hden_pos, hγpos, hγ, mul_nonneg (by linarith : (0:ℝ) ≤ (4/3)*|γ|) (by linarith : (0:ℝ) ≤ (1 + 1 / (2 * K)) - 1)]
  have hmul : ((1 + 1 / (2 * K)) + |γ| + 1 / 4) * (N : ℝ) / ((1 + 1 / (2 * K)) - 1 / 4)
      ≤ (N : ℝ) * ((4 / 3) * |γ| + 3) := by
    rw [mul_comm ((1 + 1 / (2 * K)) + |γ| + 1 / 4) (N : ℝ), mul_div_assoc]
    exact mul_le_mul_of_nonneg_left hratio hN0
  have hinner : 2 + (N : ℝ)
      + ((1 + 1 / (2 * K)) + |γ| + 1 / 4) * (N : ℝ) / ((1 + 1 / (2 * K)) - 1 / 4)
      ≤ (N : ℝ) * (2 * |γ| + 6) := by
    nlinarith [hmul, hN1, hγ, hN0, mul_nonneg hN0 hγpos.le]
  have hInnerPos : 0 < 2 + (N : ℝ)
      + ((1 + 1 / (2 * K)) + |γ| + 1 / 4) * (N : ℝ) / ((1 + 1 / (2 * K)) - 1 / 4) := by
    have hq : 0 ≤ ((1 + 1 / (2 * K)) + |γ| + 1 / 4) * (N : ℝ) / ((1 + 1 / (2 * K)) - 1 / 4) :=
      div_nonneg (mul_nonneg (by linarith) hN0) hden_pos.le
    linarith
  have hKfac : 0 < 3 + 2 * K := by linarith
  have hTargetPos : 0 < (N : ℝ) * (2 * |γ| + 6) := by
    apply mul_pos (by linarith) (by linarith)
  have hlog : Real.log ((2 + (N : ℝ)
        + ((1 + 1 / (2 * K)) + |γ| + 1 / 4) * (N : ℝ) / ((1 + 1 / (2 * K)) - 1 / 4)) * (3 + 2 * K))
      ≤ Real.log ((N : ℝ) * (2 * |γ| + 6) * (3 + 2 * K)) := by
    apply Real.log_le_log (by positivity)
    exact mul_le_mul_of_nonneg_right hinner hKfac.le
  rw [hinv]
  have hlhs : 4 * Real.log ((2 + (N : ℝ)
        + ((1 + 1 / (2 * K)) + |γ| + 1 / 4) * (N : ℝ) / ((1 + 1 / (2 * K)) - 1 / 4)) * (3 + 2 * K))
        / (1 / 4 / 2)
      = 32 * Real.log ((2 + (N : ℝ)
        + ((1 + 1 / (2 * K)) + |γ| + 1 / 4) * (N : ℝ) / ((1 + 1 / (2 * K)) - 1 / 4)) * (3 + 2 * K)) := by
    ring
  rw [hlhs]
  linarith [hlog]

#print axioms hBexp_rhs_le_LFunction
