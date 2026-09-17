import SmallTZeroFreeStrip
import Mathlib.NumberTheory.LSeries.Dirichlet

set_option maxHeartbeats 1000000

/-!
# §1617 — a zero-free disc around `s = 1` (discharges the §1616 hypotheses)

There is `R > 0` such that the closed disc `closedBall 1 R` contains no `ζ`-zero (except possibly `s=1` itself, which
is a pole) and avoids `0` and `-1`.  Take `R = min((1-σ₀)/2, 1/2)` with `σ₀` the small-|t| zero-free threshold
(§1573): every `z` in the disc has `Re z ≥ 1-R > σ₀` and `|Im z| ≤ R ≤ 1`, so `z` with `Re z ≤ 1` is covered by the
strip §1573 and `z` with `Re z > 1` by `riemannZeta_ne_zero_of_one_lt_re`; `0` and `-1` sit at distance `≥ 1 > R`.
This makes the residue-deposit §1616 unconditional for small `R`.  EXIT0 / std-axioms.
-/

open Complex

theorem exists_zero_free_ball :
    ∃ R : ℝ, 0 < R ∧ (∀ z ∈ Metric.closedBall (1 : ℂ) R, z ≠ 1 → riemannZeta z ≠ 0)
      ∧ (∀ z ∈ Metric.closedBall (1 : ℂ) R, z ≠ 0)
      ∧ (∀ z ∈ Metric.closedBall (1 : ℂ) R, z ≠ -1) := by
  obtain ⟨σ₀, hσ₀, hzf⟩ := exists_small_t_zero_free_strip
  refine ⟨min ((1 - σ₀) / 2) (1 / 2), lt_min (by linarith) (by norm_num), ?_, ?_, ?_⟩
  · -- no zeros except possibly at 1
    intro z hz hz1
    rw [Metric.mem_closedBall, Complex.dist_eq] at hz
    have hR2 : min ((1 - σ₀) / 2) (1 / 2) ≤ (1 - σ₀) / 2 := min_le_left _ _
    have hR2' : min ((1 - σ₀) / 2) (1 / 2) ≤ 1 / 2 := min_le_right _ _
    have hre : |z.re - 1| ≤ ‖z - 1‖ := by
      have := Complex.abs_re_le_norm (z - 1); rwa [Complex.sub_re, Complex.one_re] at this
    have him : |z.im| ≤ ‖z - 1‖ := by
      have := Complex.abs_im_le_norm (z - 1); rwa [Complex.sub_im, Complex.one_im, sub_zero] at this
    have hzre_lb : σ₀ < z.re := by
      have h1 : |z.re - 1| ≤ (1 - σ₀) / 2 := le_trans hre (le_trans hz hR2)
      have h2 := (abs_le.mp h1).1
      linarith
    have hzim : |z.im| ≤ 1 := by
      have : |z.im| ≤ 1 / 2 := le_trans him (le_trans hz hR2')
      linarith
    by_cases hle : z.re ≤ 1
    · exact hzf z hzre_lb hle hzim hz1
    · exact riemannZeta_ne_zero_of_one_lt_re (not_le.mp hle)
  · -- avoids 0
    intro z hz
    rw [Metric.mem_closedBall, Complex.dist_eq] at hz
    intro h; rw [h] at hz
    simp only [zero_sub, norm_neg, norm_one] at hz
    have := min_le_right ((1 - σ₀) / 2) (1 / 2)
    linarith
  · -- avoids -1
    intro z hz
    rw [Metric.mem_closedBall, Complex.dist_eq] at hz
    intro h; rw [h] at hz
    have hnorm : ‖(-1 : ℂ) - 1‖ = 2 := by
      norm_num [show ((-1 : ℂ) - 1) = (-2 : ℂ) by ring]
    rw [hnorm] at hz
    have := min_le_right ((1 - σ₀) / 2) (1 / 2)
    linarith

#print axioms exists_zero_free_ball
