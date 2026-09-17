import ExistsBranchLogBall
import HasDerivBranch

/-!
# §1324 — Holomorphic log-branch on the ball: `∃ f, HasDerivAt f (logDeriv h) & exp∘f = h`

Assembling §1322 (continuous branch on the ball) with §1323 (continuous branch of a holomorphic function is
holomorphic, deriv `= logDeriv`).  For `h` differentiable and nonvanishing on `ball 0 R`, there is a branch `f`
of `log h` that is holomorphic at every point of the ball with `HasDerivAt f (logDeriv h z) z` and `exp∘f = h`.
This packages the branch construction into exactly the object Borel–Carathéodory §1318 consumes (via §1320's
`Re f = log‖h‖` and this derivative).  EXIT0 / std-axioms.
-/

open Complex Metric Filter Topology

theorem exists_holo_branch_on_ball {R : ℝ} (hR : 0 < R) {h : ℂ → ℂ}
    (hh : DifferentiableOn ℂ h (ball 0 R)) (hh0 : ∀ z ∈ ball 0 R, h z ≠ 0) :
    ∃ f : ℂ → ℂ, (∀ z ∈ ball (0 : ℂ) R, HasDerivAt f (logDeriv h z) z) ∧
      (∀ z ∈ ball (0 : ℂ) R, Complex.exp (f z) = h z) := by
  obtain ⟨f, hfc, hfe⟩ := exists_branch_log_on_ball hR hh.continuousOn hh0
  refine ⟨f, fun z hz => ?_, hfe⟩
  have hnhds : ball (0 : ℂ) R ∈ 𝓝 z := isOpen_ball.mem_nhds hz
  have hφc : ContinuousAt f z := (hfc z hz).continuousAt hnhds
  have hψd : HasDerivAt h (deriv h z) z := ((hh z hz).differentiableAt hnhds).hasDerivAt
  have hcomp : ∀ᶠ w in 𝓝 z, Complex.exp (f w) = h w :=
    Filter.eventually_of_mem hnhds (fun w hw => hfe w hw)
  have hd := hasDerivAt_of_exp_comp hφc hψd hcomp
  rw [logDeriv_apply]
  exact hd

#print axioms exists_holo_branch_on_ball
