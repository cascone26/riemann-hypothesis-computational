import NormLogDerivBall

set_option maxHeartbeats 1000000

open Complex Metric

/-!
# §1670 — norm-Landau core, arbitrary center: translation of §1669 to `ball c R`

The connector points live near `1 + iT`, not the origin — this brick moves §1669's bound to any
center by the substitution `k(z) = h(z + c)` (`logDeriv k z = logDeriv h (z+c)` via
`deriv_comp_add_const`).
-/

theorem norm_logDeriv_le_of_forall_norm_le_center {c : ℂ} {R M' : ℝ} {h : ℂ → ℂ}
    (hR : 0 < R) (hM' : 0 < M')
    (hh : DifferentiableOn ℂ h (ball c R)) (hh0 : ∀ z ∈ ball c R, h z ≠ 0)
    (hgrow : ∀ z ∈ ball c R, ‖h z‖ ≤ Real.exp M' * ‖h c‖) :
    ∀ s ∈ closedBall c (R / 4), ‖logDeriv h s‖ ≤ 8 * M' / R := by
  set k : ℂ → ℂ := fun z => h (z + c) with hkdef
  have hmap : ∀ z ∈ ball (0 : ℂ) R, z + c ∈ ball c R := by
    intro z hz
    rw [mem_ball, dist_zero_right] at hz
    rw [mem_ball, dist_eq_norm]
    simpa using hz
  have hkd : DifferentiableOn ℂ k (ball 0 R) := by
    intro z hz
    exact ((hh (z + c) (hmap z hz)).differentiableAt
      (isOpen_ball.mem_nhds (hmap z hz))).comp z
      ((differentiableAt_id.add_const c)) |>.differentiableWithinAt
  have hk0 : ∀ z ∈ ball (0 : ℂ) R, k z ≠ 0 := fun z hz => hh0 (z + c) (hmap z hz)
  have hkc : k 0 = h c := by rw [hkdef]; simp
  have hkgrow : ∀ z ∈ ball (0 : ℂ) R, ‖k z‖ ≤ Real.exp M' * ‖k 0‖ := by
    intro z hz
    rw [hkc]
    exact hgrow (z + c) (hmap z hz)
  have hklog : ∀ z : ℂ, logDeriv k z = logDeriv h (z + c) := by
    intro z
    rw [logDeriv_apply, logDeriv_apply, hkdef]
    congr 1
    exact deriv_comp_add_const h c z
  intro s hs
  have hs' : s - c ∈ closedBall (0 : ℂ) (R / 4) := by
    rw [mem_closedBall, dist_zero_right]
    rw [mem_closedBall, dist_eq_norm] at hs
    exact hs
  have := norm_logDeriv_le_of_forall_norm_le hR hM' hkd hk0 hkgrow (s - c) hs'
  rwa [hklog, sub_add_cancel] at this

#print axioms norm_logDeriv_le_of_forall_norm_le_center
