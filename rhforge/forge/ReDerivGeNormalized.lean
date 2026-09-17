import ReDerivGeOfReLe

/-!
# §1367 — normalized Borel–Carathéodory derivative bound (kills the `‖f 0‖` term)

§1318 (`re_deriv_ge_of_re_le`) carries an uncontrolled `6‖f 0‖` term.  Applying it to the *recentred* function
`f - f 0` (which vanishes at `0`, so its `‖·0‖` term drops) gives the standard, clean Borel–Carathéodory form:

  `-(4 (M - Re f 0) / R) ≤ Re (deriv f 0)`,

i.e. the derivative bound depends only on `M - Re f 0 = sup Re f - Re f(0)`.  Since `deriv` and `logDeriv` are
invariant under the additive/constant normalisation, this is the version that discharges §1355's `hgB` with an
absolute `B` once `sup‖g‖/‖g c‖` on the disk is bounded (via Jensen §1361–§1366).  EXIT0 / std-axioms.
-/

open Complex

theorem re_deriv_ge_of_re_le_normalized {f : ℂ → ℂ} {R M : ℝ} (hR : 0 < R)
    (hM' : 0 < M - (f 0).re)
    (hf : DifferentiableOn ℂ f (Metric.ball 0 R))
    (hre : ∀ z ∈ Metric.ball (0 : ℂ) R, (f z).re ≤ M) :
    -(4 * (M - (f 0).re) / R) ≤ (deriv f 0).re := by
  set g : ℂ → ℂ := fun z => f z - f 0 with hg
  have hgd : DifferentiableOn ℂ g (Metric.ball 0 R) := hf.sub_const _
  have hgre : ∀ z ∈ Metric.ball (0 : ℂ) R, (g z).re ≤ M - (f 0).re := by
    intro z hz
    simp only [hg, Complex.sub_re]
    linarith [hre z hz]
  have hg0 : g 0 = 0 := by simp [hg]
  have hd := re_deriv_ge_of_re_le hR hM' hgd hgre
  rw [hg0, norm_zero, mul_zero, add_zero] at hd
  have hderiv : deriv g 0 = deriv f 0 := by
    rw [hg]; exact deriv_sub_const _
  rwa [hderiv] at hd

#print axioms re_deriv_ge_of_re_le_normalized
