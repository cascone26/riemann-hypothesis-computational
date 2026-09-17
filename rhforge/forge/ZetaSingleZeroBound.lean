import ZetaZeroCountBound
import AnalyticAtZeta

/-!
# §1364 — single-zero bound for ζ (extracting one term from §1363)

On a pole-free ball (`1 ∉ closedBall c R`) every summand of the zero-count sum §1363 is `≥ 0` (`divisor ≥ 0`
since ζ is analytic there; `log(R/‖c-u‖) ≥ 0` since support points lie in the ball).  Hence by
`single_le_finsum` the contribution of any single point `ρ` is bounded by the whole sum, which §1363 bounds:

  `divisor ζ ρ · log(R/‖c-ρ‖) ≤ log M - divisor ζ c · log R - log‖trailingCoeff ζ c‖`.

For a zero `ρ` (`divisor ρ ≥ 1`) near the centre this forces `log(R/‖c-ρ‖) ≍ log M`, i.e. zeros keep their
distance — the quantitative input to `‖ζ/P‖ ≤ M'` and §1355's `hgB`.  EXIT0 / std-axioms.
-/

open Complex Metric Real

theorem riemannZeta_single_zero_le {c ρ : ℂ} {R M : ℝ} (hR : 0 < R) (hM : 1 ≤ M)
    (hball : (1 : ℂ) ∉ closedBall c |R|)
    (hbound : ∀ s ∈ sphere c |R|, ‖riemannZeta s‖ ≤ M) :
    (MeromorphicOn.divisor riemannZeta (closedBall c |R|) ρ : ℝ) * Real.log (R * ‖c - ρ‖⁻¹)
      ≤ Real.log M
        - (MeromorphicOn.divisor riemannZeta (closedBall c |R|) c : ℝ) * Real.log R
        - Real.log ‖meromorphicTrailingCoeffAt riemannZeta c‖ := by
  set D := MeromorphicOn.divisor riemannZeta (closedBall c |R|) with hD
  set f : ℂ → ℝ := fun u => (D u : ℝ) * Real.log (R * ‖c - u‖⁻¹) with hf
  have hmero : MeromorphicOn riemannZeta (closedBall c |R|) :=
    meromorphicOn_riemannZeta.mono_set (Set.subset_univ _)
  have hnn : ∀ u, 0 ≤ f u := by
    intro u
    rcases eq_or_ne (D u) 0 with hd | hd
    · simp [hf, hd]
    · have humem : u ∈ closedBall c |R| := D.supportWithinDomain hd
      have hu1 : u ≠ 1 := fun h => hball (h ▸ humem)
      have hdnn : 0 ≤ D u := by
        have hord := (analyticAt_riemannZeta hu1).meromorphicOrderAt_nonneg
        rw [hD, MeromorphicOn.divisor_apply hmero humem]
        rcases eq_or_ne (meromorphicOrderAt riemannZeta u) ⊤ with h | h
        · rw [h]; simp
        · lift meromorphicOrderAt riemannZeta u to ℤ using h with n; simpa using hord
      have hlog : 0 ≤ Real.log (R * ‖c - u‖⁻¹) := by
        rcases eq_or_ne u c with rfl | huc
        · simp
        · have hcu : 0 < ‖c - u‖ := by
            rw [norm_pos_iff, sub_ne_zero]; exact fun h => huc h.symm
          have hle : ‖c - u‖ ≤ R := by
            rw [mem_closedBall, Complex.dist_eq] at humem
            rw [show ‖c - u‖ = ‖u - c‖ from by rw [norm_sub_rev]]
            calc ‖u - c‖ ≤ |R| := humem
              _ = R := abs_of_pos hR
          apply Real.log_nonneg
          rw [← div_eq_mul_inv]
          exact (one_le_div hcu).mpr hle
      exact mul_nonneg (by exact_mod_cast hdnn) hlog
  have hDfin : (Function.support (fun u => D u)).Finite := D.finiteSupport (isCompact_closedBall c |R|)
  have hsub : Function.support f ⊆ Function.support (fun u => D u) := by
    intro u hu
    simp only [Function.mem_support] at hu ⊢
    intro h; apply hu; simp [hf, h]
  have hfin : Function.HasFiniteSupport f := hDfin.subset hsub
  have hsingle := single_le_finsum ρ hfin hnn
  exact le_trans hsingle (riemannZeta_zero_count_le hR.ne' hM hbound)

#print axioms riemannZeta_single_zero_le
