import ReLogDerivGeNormalized

/-!
# §1369 — center-`c` normalized Borel–Carathéodory bound for `logDeriv`

The center-`c` translate of §1368 (as §1352 was for §1325), via the shift `sh w = c + w`.  For `h` differentiable
and nonvanishing on `ball c R` with `‖h‖ ≤ M'` and `‖h c‖ < M'`:

  `-(4 (log M' - log‖h c‖) / R) ≤ Re (logDeriv h c)`.

Evaluated at the actual DLVP test centre `c = σ + iγ`, and with `h = ζ/P` the nonvanishing factor, this is
exactly §1355's `hgB` once `M'` (the disk sup of `‖ζ/P‖`) is supplied by Jensen (§1361–§1366).  No `‖f 0‖` term.
EXIT0 / std-axioms.
-/

open Complex Metric Real

theorem re_logDeriv_ge_of_bounded_normalized_at {c : ℂ} {R M' : ℝ} (hR : 0 < R) {h : ℂ → ℂ}
    (hh : DifferentiableOn ℂ h (ball c R)) (hh0 : ∀ z ∈ ball c R, h z ≠ 0)
    (hb : ∀ z ∈ ball c R, ‖h z‖ ≤ M') (hh0lt : ‖h c‖ < M') :
    -(4 * (Real.log M' - Real.log ‖h c‖) / R) ≤ (logDeriv h c).re := by
  set sh : ℂ → ℂ := fun w => c + w with hsh
  have hmaps : Set.MapsTo sh (ball 0 R) (ball c R) := by
    intro w hw
    rw [mem_ball, dist_zero_right] at hw
    rw [mem_ball, Complex.dist_eq]
    simpa [hsh] using hw
  have hshdiff : ∀ x, DifferentiableAt ℂ sh x := by
    intro x; rw [hsh]; fun_prop
  have hshon : DifferentiableOn ℂ sh (ball (0 : ℂ) R) :=
    fun x _ => (hshdiff x).differentiableWithinAt
  have hh' : DifferentiableOn ℂ (h ∘ sh) (ball 0 R) := hh.comp hshon hmaps
  have hh0' : ∀ z ∈ ball (0 : ℂ) R, (h ∘ sh) z ≠ 0 := fun z hz => hh0 _ (hmaps hz)
  have hb' : ∀ z ∈ ball (0 : ℂ) R, ‖(h ∘ sh) z‖ ≤ M' := fun z hz => hb _ (hmaps hz)
  have h00 : (h ∘ sh) 0 = h c := by simp [hsh]
  have hb0lt : ‖(h ∘ sh) 0‖ < M' := by rw [h00]; exact hh0lt
  have hdc : DifferentiableAt ℂ h c :=
    (hh c (mem_ball_self hR)).differentiableAt (isOpen_ball.mem_nhds (mem_ball_self hR))
  have hld : logDeriv (h ∘ sh) 0 = logDeriv h c := by
    have hcomp := logDeriv_comp (by simpa [hsh] using hdc) (hshdiff 0)
    have hsh0 : sh 0 = c := by simp [hsh]
    have hderiv : deriv sh 0 = 1 := by
      simpa [hsh] using ((hasDerivAt_id (0 : ℂ)).const_add c).deriv
    rw [hcomp, hsh0, hderiv, mul_one]
  have hbc := re_logDeriv_ge_of_bounded_normalized hR hh' hh0' hb' hb0lt
  rw [h00, hld] at hbc
  exact hbc

#print axioms re_logDeriv_ge_of_bounded_normalized_at
