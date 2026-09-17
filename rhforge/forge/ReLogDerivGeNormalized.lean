import ExistsHoloBranchBall
import ReDerivGeNormalized
import ReEqLogNormOfExp

/-!
# §1368 — normalized Borel–Carathéodory bound for `logDeriv` (the clean `hgB` tool)

The normalized companion of §1325: threading the holomorphic branch (§1324) through the *normalized* BC bound
§1367 instead of §1318 removes the uncontrolled `‖f 0‖` term.  For `h` differentiable and nonvanishing on
`ball 0 R` with `‖h‖ ≤ M'` and `‖h 0‖ < M'`:

  `-(4 (log M' - log‖h 0‖) / R) ≤ Re (logDeriv h 0)`.

The bound depends only on `log M' - log‖h 0‖ = log(M'/‖h 0‖)` — an absolute `≍ log|t|` quantity once the ζ/P
sup-bound and centre value are supplied by Jensen (§1361–§1366).  This is exactly the shape §1355's `hgB` needs.
EXIT0 / std-axioms.
-/

open Complex Metric Real

theorem re_logDeriv_ge_of_bounded_normalized {R M' : ℝ} (hR : 0 < R) {h : ℂ → ℂ}
    (hh : DifferentiableOn ℂ h (ball 0 R)) (hh0 : ∀ z ∈ ball 0 R, h z ≠ 0)
    (hb : ∀ z ∈ ball 0 R, ‖h z‖ ≤ M') (hh0lt : ‖h 0‖ < M') :
    -(4 * (Real.log M' - Real.log ‖h 0‖) / R) ≤ (logDeriv h 0).re := by
  obtain ⟨f, hfd, hfe⟩ := exists_holo_branch_on_ball hR hh hh0
  have h0ball : (0 : ℂ) ∈ ball 0 R := mem_ball_self hR
  have hh00 : h 0 ≠ 0 := hh0 0 h0ball
  have hfdiff : DifferentiableOn ℂ f (ball 0 R) :=
    fun z hz => (hfd z hz).differentiableAt.differentiableWithinAt
  have hf0 : (f 0).re = Real.log ‖h 0‖ := re_eq_log_norm_of_exp (hfe 0 h0ball)
  have hre : ∀ z ∈ ball (0 : ℂ) R, (f z).re ≤ Real.log M' := by
    intro z hz
    rw [re_eq_log_norm_of_exp (hfe z hz)]
    exact Real.log_le_log (norm_pos_iff.mpr (hh0 z hz)) (hb z hz)
  have hM0' : 0 < Real.log M' - (f 0).re := by
    rw [hf0]
    have := Real.log_lt_log (norm_pos_iff.mpr hh00) hh0lt
    linarith
  have hbc := re_deriv_ge_of_re_le_normalized hR hM0' hfdiff hre
  rw [hf0] at hbc
  rwa [(hfd 0 h0ball).deriv] at hbc

#print axioms re_logDeriv_ge_of_bounded_normalized
