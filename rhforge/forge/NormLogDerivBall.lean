import ExistsHoloBranchBall
import Mathlib.Analysis.Complex.BorelCaratheodory
import Mathlib.Analysis.Complex.Liouville

set_option maxHeartbeats 1000000

open Complex Metric

/-!
# §1669 ★★ — NORM-LANDAU CORE (zero-free case): `‖logDeriv h‖ ≤ 8M'/R` on the quarter-ball

The tower's Landau decomposition (§ HdecompHgB) controls only `Re g ≥ -B` — all DLVP needed. The
effective-PNT connector bound needs the NORM of the zero-free factor's log-derivative. This brick
supplies it: `h` holomorphic and nonvanishing on `ball 0 R` with the growth bound
`‖h z‖ ≤ e^{M'}·‖h 0‖` has `‖logDeriv h s‖ ≤ 8M'/R` for `‖s‖ ≤ R/4`. Chain, fully in-library:
our §1323 `exists_holo_branch_on_ball` (analytic branch `f`, `f' = logDeriv h`, `exp f = h`) →
`Re(f - f 0) = log(‖h‖/‖h 0‖) ≤ M'` → Mathlib `borelCaratheodory_zero` (‖f-f0‖ ≤ 2M' on the
half-ball) → Mathlib Cauchy estimate `norm_deriv_le_of_forall_mem_sphere_norm_le` on spheres of
radius `R/4` centered in the quarter-ball.
-/

theorem norm_logDeriv_le_of_forall_norm_le {R M' : ℝ} {h : ℂ → ℂ} (hR : 0 < R) (hM' : 0 < M')
    (hh : DifferentiableOn ℂ h (ball 0 R)) (hh0 : ∀ z ∈ ball 0 R, h z ≠ 0)
    (hgrow : ∀ z ∈ ball 0 R, ‖h z‖ ≤ Real.exp M' * ‖h 0‖) :
    ∀ s ∈ closedBall (0 : ℂ) (R / 4), ‖logDeriv h s‖ ≤ 8 * M' / R := by
  have h0mem : (0 : ℂ) ∈ ball (0 : ℂ) R := mem_ball_self hR
  have h0pos : 0 < ‖h 0‖ := norm_pos_iff.mpr (hh0 0 h0mem)
  obtain ⟨f, hfd, hfe⟩ := exists_holo_branch_on_ball hR hh hh0
  set F : ℂ → ℂ := fun z => f z - f 0 with hFdef
  -- Re f z = log ‖h z‖ on the ball
  have href : ∀ z ∈ ball (0 : ℂ) R, (f z).re = Real.log ‖h z‖ := by
    intro z hz
    have : ‖h z‖ = Real.exp (f z).re := by
      rw [← hfe z hz, Complex.norm_exp]
    rw [this, Real.log_exp]
  have hRe : ∀ z ∈ ball (0 : ℂ) R, (F z).re ≤ M' := by
    intro z hz
    have hz' : (0 : ℝ) < ‖h z‖ := norm_pos_iff.mpr (hh0 z hz)
    have hlog : Real.log ‖h z‖ ≤ M' + Real.log ‖h 0‖ := by
      have := Real.log_le_log hz' (hgrow z hz)
      rwa [Real.log_mul (Real.exp_ne_zero _) h0pos.ne', Real.log_exp] at this
    show (f z - f 0).re ≤ M'
    rw [Complex.sub_re, href z hz, href 0 h0mem]
    linarith
  have hFdiff : DifferentiableOn ℂ F (ball 0 R) := by
    intro z hz
    exact (((hfd z hz).differentiableAt).sub_const (f 0)).differentiableWithinAt
  -- Borel-Carathéodory: ‖F‖ ≤ 2M' on the closed half-ball
  have hBC : ∀ z ∈ closedBall (0 : ℂ) (R / 2), ‖F z‖ ≤ 2 * M' := by
    intro z hz
    rw [mem_closedBall, dist_zero_right] at hz
    have hzball : z ∈ ball (0 : ℂ) R := by
      rw [mem_ball, dist_zero_right]
      linarith
    have hbc := Complex.borelCaratheodory_zero hM' hFdiff
      (fun u hu => hRe u hu) hR hzball (by simp [hFdef])
    have hdenom : 0 < R - ‖z‖ := by
      rw [mem_ball, dist_zero_right] at hzball
      linarith
    calc ‖F z‖ ≤ 2 * M' * ‖z‖ / (R - ‖z‖) := hbc
      _ ≤ 2 * M' := by
          rw [div_le_iff₀ hdenom]
          nlinarith [norm_nonneg z]
  -- Cauchy estimate on spheres of radius R/4 around points of the quarter-ball
  intro s hs
  rw [mem_closedBall, dist_zero_right] at hs
  have hsub : closedBall s (R / 4) ⊆ closedBall (0 : ℂ) (R / 2) := by
    intro z hz
    rw [mem_closedBall, dist_zero_right]
    rw [mem_closedBall, dist_eq_norm] at hz
    calc ‖z‖ = ‖z - s + s‖ := by ring_nf
      _ ≤ ‖z - s‖ + ‖s‖ := norm_add_le _ _
      _ ≤ R / 4 + R / 4 := by linarith
      _ = R / 2 := by ring
  have hsubR : closedBall s (R / 4) ⊆ ball (0 : ℂ) R := by
    intro z hz
    have := hsub hz
    rw [mem_closedBall, dist_zero_right] at this
    rw [mem_ball, dist_zero_right]
    linarith
  have hFcl : DiffContOnCl ℂ F (ball s (R / 4)) := by
    have hd : DifferentiableOn ℂ F (closedBall s (R / 4)) := hFdiff.mono hsubR
    rw [← closure_ball s (by positivity : (R / 4 : ℝ) ≠ 0)] at hd
    exact hd.diffContOnCl
  have hsphere : ∀ z ∈ sphere s (R / 4), ‖F z‖ ≤ 2 * M' := by
    intro z hz
    exact hBC z (hsub (sphere_subset_closedBall hz))
  have hcauchy := Complex.norm_deriv_le_of_forall_mem_sphere_norm_le
    (by positivity : (0 : ℝ) < R / 4) hFcl hsphere
  -- deriv F s = logDeriv h s
  have hsball : s ∈ ball (0 : ℂ) R := by
    rw [mem_ball, dist_zero_right]
    linarith
  have hderivF : deriv F s = logDeriv h s := by
    have h1 : HasDerivAt F (logDeriv h s) s := by
      simpa [hFdef] using (hfd s hsball).sub_const (f 0)
    exact h1.deriv
  rw [hderivF] at hcauchy
  calc ‖logDeriv h s‖ ≤ 2 * M' / (R / 4) := hcauchy
    _ = 8 * M' / R := by
        field_simp
        ring

#print axioms norm_logDeriv_le_of_forall_norm_le
