import ReLogDerivGe

/-!
# §1352 — center-`c` Borel–Carathéodory logDeriv bound (translate of §1325)

§1325 (`re_logDeriv_ge_of_bounded`) bounds `Re (logDeriv h 0)` on `ball 0 R` — at the *center* only.  The grand
assembly (§1336's `hgB`) needs the bound at the actual DLVP test point `s₀ = σ + iγ`, which is the *center* of the
Borel–Carathéodory disk but is not `0`.  Translating by `sh w = c + w` (`h ∘ sh` on `ball 0 R`, `MapsTo` onto
`ball c R`) and using `logDeriv_comp` (chain rule: `logDeriv (h ∘ sh) 0 = logDeriv h (sh 0) · deriv sh 0`, and
`sh 0 = c`, `deriv sh 0 = 1`) carries §1325 to an arbitrary center `c`.  So for `h` differentiable, nonvanishing,
`‖h‖ ≤ M'` on `ball c R`:  `-((4 log M' + 6‖f 0‖)/R) ≤ (logDeriv h c).re` for a branch `f` with `exp (f 0) = h c`.
Directly consumable at `c = σ + iγ`.  EXIT0 / std-axioms.
-/

open Complex Metric

theorem re_logDeriv_ge_of_bounded_at {c : ℂ} {R M' : ℝ} (hR : 0 < R) (hM' : 1 < M')
    {h : ℂ → ℂ} (hh : DifferentiableOn ℂ h (ball c R)) (hh0 : ∀ z ∈ ball c R, h z ≠ 0)
    (hb : ∀ z ∈ ball c R, ‖h z‖ ≤ M') :
    ∃ f : ℂ → ℂ, Complex.exp (f 0) = h c ∧
      -((4 * Real.log M' + 6 * ‖f 0‖) / R) ≤ (logDeriv h c).re := by
  set sh : ℂ → ℂ := fun w => c + w with hsh
  have hmaps : Set.MapsTo sh (ball 0 R) (ball c R) := by
    intro w hw
    rw [mem_ball, dist_zero_right] at hw
    rw [mem_ball, dist_eq]
    simpa [hsh] using hw
  have hshdiff : ∀ x, DifferentiableAt ℂ sh x := by
    intro x; rw [hsh]; fun_prop
  have hshon : DifferentiableOn ℂ sh (ball (0 : ℂ) R) :=
    fun x _ => (hshdiff x).differentiableWithinAt
  have hh' : DifferentiableOn ℂ (h ∘ sh) (ball 0 R) := hh.comp hshon hmaps
  have hh0' : ∀ z ∈ ball (0 : ℂ) R, (h ∘ sh) z ≠ 0 := fun z hz => hh0 _ (hmaps hz)
  have hb' : ∀ z ∈ ball (0 : ℂ) R, ‖(h ∘ sh) z‖ ≤ M' := fun z hz => hb _ (hmaps hz)
  obtain ⟨f, hfe, hbc⟩ := re_logDeriv_ge_of_bounded hR hM' hh' hh0' hb'
  have h00 : (h ∘ sh) 0 = h c := by simp [hsh]
  have hdc : DifferentiableAt ℂ h c :=
    (hh c (mem_ball_self hR)).differentiableAt (isOpen_ball.mem_nhds (mem_ball_self hR))
  have hld : logDeriv (h ∘ sh) 0 = logDeriv h c := by
    have hcomp := logDeriv_comp (by simpa [hsh] using hdc) (hshdiff 0)
    have hsh0 : sh 0 = c := by simp [hsh]
    have hderiv : deriv sh 0 = 1 := by
      simpa [hsh] using ((hasDerivAt_id (0 : ℂ)).const_add c).deriv
    rw [hcomp, hsh0, hderiv, mul_one]
  exact ⟨f, hfe.trans h00, by rw [hld] at hbc; exact hbc⟩

#print axioms re_logDeriv_ge_of_bounded_at
