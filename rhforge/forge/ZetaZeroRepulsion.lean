import ZetaSingleZeroBound

/-!
# §1366 — zero-repulsion bound for ζ (multiplicity-free form of §1364)

A genuine zero `ρ` (`divisor ρ ≥ 1`) satisfies `log(R/‖c-ρ‖) ≤ divisor ρ · log(R/‖c-ρ‖)` (since the log factor
is `≥ 0`), so §1364 gives directly

  `log(R/‖c-ρ‖) ≤ log M - divisor ζ c · log R - log‖trailingCoeff ζ c‖`.

Equivalently `R/‖c-ρ‖ ≤ M · …`, i.e. `ρ` cannot approach the centre faster than the (`≍ log|t|`) bound allows —
the clean zero-repulsion statement feeding the factor bound `‖ζ/P‖ ≤ M'`.  EXIT0 / std-axioms.
-/

open Complex Metric Real

theorem riemannZeta_zero_log_dist_le {c ρ : ℂ} {R M : ℝ} (hR : 0 < R) (hM : 1 ≤ M)
    (hball : (1 : ℂ) ∉ closedBall c |R|)
    (hbound : ∀ s ∈ sphere c |R|, ‖riemannZeta s‖ ≤ M)
    (hρmem : ρ ∈ closedBall c |R|) (hρc : ρ ≠ c)
    (hdiv1 : 1 ≤ MeromorphicOn.divisor riemannZeta (closedBall c |R|) ρ) :
    Real.log (R * ‖c - ρ‖⁻¹)
      ≤ Real.log M
        - (MeromorphicOn.divisor riemannZeta (closedBall c |R|) c : ℝ) * Real.log R
        - Real.log ‖meromorphicTrailingCoeffAt riemannZeta c‖ := by
  have hcu : 0 < ‖c - ρ‖ := by
    rw [norm_pos_iff, sub_ne_zero]; exact fun h => hρc h.symm
  have hle : ‖c - ρ‖ ≤ R := by
    rw [mem_closedBall, Complex.dist_eq] at hρmem
    rw [show ‖c - ρ‖ = ‖ρ - c‖ from by rw [norm_sub_rev]]
    calc ‖ρ - c‖ ≤ |R| := hρmem
      _ = R := abs_of_pos hR
  have hlog : 0 ≤ Real.log (R * ‖c - ρ‖⁻¹) := by
    apply Real.log_nonneg; rw [← div_eq_mul_inv]; exact (one_le_div hcu).mpr hle
  have hd1 : (1 : ℝ) ≤ (MeromorphicOn.divisor riemannZeta (closedBall c |R|) ρ : ℝ) := by
    exact_mod_cast hdiv1
  calc Real.log (R * ‖c - ρ‖⁻¹)
      ≤ (MeromorphicOn.divisor riemannZeta (closedBall c |R|) ρ : ℝ) * Real.log (R * ‖c - ρ‖⁻¹) :=
        le_mul_of_one_le_left hlog hd1
    _ ≤ _ := riemannZeta_single_zero_le hR hM hball hbound

#print axioms riemannZeta_zero_log_dist_le
