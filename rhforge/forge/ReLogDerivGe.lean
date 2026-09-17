import ExistsHoloBranchBall
import ReDerivGeOfReLe
import ReEqLogNormOfExp

/-!
# §1325 ★ Borel–Carathéodory bound on `Re(logDeriv h 0)` for bounded nonvanishing `h`

The full Borel–Carathéodory estimate, assembled: for `h` differentiable and nonvanishing on `ball 0 R` with
`‖h‖ ≤ M'` (`M' > 1`), take the holomorphic branch `f` of `log h` (§1324); then `Re f = log‖h‖ ≤ log M'` (§1320),
`f` differentiable, so §1318 gives `Re(deriv f 0) ≥ -(4·log M' + 6‖f 0‖)/R`, and `deriv f 0 = logDeriv h 0`:

  `Re(logDeriv h 0) ≥ -(4·log M' + 6‖f 0‖)/R`.

Applied to `h = ζ/P` on the disk (bounded by the ball bound §1130, `log M' ≍ log|t|`), this is exactly the
`Re g ≥ -B` input (`g = logDeriv(ζ/P)`) that §1314 turns into `hzero` and §1311 into the region.  EXIT0 / std-axioms.
-/

open Complex Metric

theorem re_logDeriv_ge_of_bounded {R M' : ℝ} (hR : 0 < R) (hM' : 1 < M') {h : ℂ → ℂ}
    (hh : DifferentiableOn ℂ h (ball 0 R)) (hh0 : ∀ z ∈ ball 0 R, h z ≠ 0)
    (hb : ∀ z ∈ ball 0 R, ‖h z‖ ≤ M') :
    ∃ f : ℂ → ℂ, Complex.exp (f 0) = h 0 ∧
      -((4 * Real.log M' + 6 * ‖f 0‖) / R) ≤ (logDeriv h 0).re := by
  obtain ⟨f, hfd, hfe⟩ := exists_holo_branch_on_ball hR hh hh0
  have h0ball : (0 : ℂ) ∈ ball 0 R := mem_ball_self hR
  refine ⟨f, hfe 0 h0ball, ?_⟩
  have hMpos : 0 < Real.log M' := Real.log_pos hM'
  have hfdiff : DifferentiableOn ℂ f (ball 0 R) :=
    fun z hz => (hfd z hz).differentiableAt.differentiableWithinAt
  have hre : ∀ z ∈ ball (0 : ℂ) R, (f z).re ≤ Real.log M' := by
    intro z hz
    rw [re_eq_log_norm_of_exp (hfe z hz)]
    exact Real.log_le_log (norm_pos_iff.mpr (hh0 z hz)) (hb z hz)
  have hbc := re_deriv_ge_of_re_le hR hMpos hfdiff hre
  rwa [(hfd 0 h0ball).deriv] at hbc

#print axioms re_logDeriv_ge_of_bounded
