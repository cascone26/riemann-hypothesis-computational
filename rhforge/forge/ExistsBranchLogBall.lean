import Mathlib.Analysis.Complex.BranchLogRoot
import Mathlib.Analysis.Normed.Module.Connected

/-!
# §1322 — Continuous branch of `log g` on a ball: existence via Mathlib's simply-connected branch

The disk `ball 0 R` is convex, hence contractible (`contractibleSpace_ball`), hence simply connected
(`SimplyConnectedSpace.ofContractible`).  So Mathlib's `exists_continuousOn_eqOn_exp_comp` applies: for any
continuous nonvanishing `g` on the ball there is a CONTINUOUS branch `f` of `log∘g` (`exp∘f = g`).  This is the
existence half of the sole remaining ζ-analytic gap (the branch of `log(ζ/P)`) — Mathlib DOES have it after all,
for continuity.  The remaining upgrade is continuous→holomorphic (exp is a local biholomorphism).  EXIT0 / std-axioms.
-/

open Complex Metric

theorem exists_branch_log_on_ball {R : ℝ} (hR : 0 < R) {g : ℂ → ℂ}
    (hgc : ContinuousOn g (ball (0 : ℂ) R)) (hg0 : ∀ z ∈ ball (0 : ℂ) R, g z ≠ 0) :
    ∃ f : ℂ → ℂ, ContinuousOn f (ball (0 : ℂ) R) ∧
      ∀ z ∈ ball (0 : ℂ) R, Complex.exp (f z) = g z := by
  haveI : ContractibleSpace (ball (0 : ℂ) R) := contractibleSpace_ball hR
  have hsc : IsSimplyConnected (ball (0 : ℂ) R) := inferInstanceAs (SimplyConnectedSpace _)
  have hU0 : (0 : ℂ) ∉ g '' (ball (0 : ℂ) R) := by
    rintro ⟨z, hz, hgz⟩; exact hg0 z hz hgz
  obtain ⟨f, hfc, hfe⟩ := Complex.exists_continuousOn_eqOn_exp_comp hsc isOpen_ball hgc hU0
  exact ⟨f, hfc, fun z hz => hfe hz⟩

#print axioms exists_branch_log_on_ball
