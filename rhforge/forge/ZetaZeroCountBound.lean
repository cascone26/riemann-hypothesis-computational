import ZetaJensen
import ZetaCircleAvgBound

/-!
# §1363 — zero-count bound for ζ (Jensen §1361 + circle bound §1362)

Combining Jensen's formula §1361 with the circle bound §1362 gives an upper bound on the divisor-weighted
zero-count sum:

  `∑ᶠ u, divisor ζ u · log(R/‖c-u‖) ≤ log M - divisor ζ c · log R - log‖trailingCoeff ζ c‖`.

On a pole-free ball every summand is `≥ 0` (`divisor ≥ 0`, `log(R/‖c-u‖) ≥ 0` for `u` in the ball), so this
bounds the contribution of any individual zero — the quantitative "zeros cannot cluster near the centre" fact
that yields `‖ζ/P‖ ≤ M'` and §1355's `hgB`.  EXIT0 / std-axioms.
-/

open Complex Metric Real

theorem riemannZeta_zero_count_le {c : ℂ} {R M : ℝ} (hR : R ≠ 0) (hM : 1 ≤ M)
    (hbound : ∀ s ∈ sphere c |R|, ‖riemannZeta s‖ ≤ M) :
    (∑ᶠ u, (MeromorphicOn.divisor riemannZeta (closedBall c |R|) u : ℝ) * Real.log (R * ‖c - u‖⁻¹))
      ≤ Real.log M
        - (MeromorphicOn.divisor riemannZeta (closedBall c |R|) c : ℝ) * Real.log R
        - Real.log ‖meromorphicTrailingCoeffAt riemannZeta c‖ := by
  have h := riemannZeta_circleAverage_log_norm (c := c) (R := R) hR
  have hle := circleAverage_log_norm_riemannZeta_le hM hbound
  rw [h] at hle
  linarith

#print axioms riemannZeta_zero_count_le
