import ZetaEStripGrowth
import ZetaThreeLines

/-!
# S2999 — hB: E·G_T is bounded on the closed critical strip (lane 1: last three-lines hypothesis)

F_T(s) = E(s)·exp((s − iT)²), E = (s−1)ζ(s).
* `norm_F_far` — 0 ≤ σ ≤ 1, |t| ≥ 1 ⇒ ‖F_T(s)‖ ≤ 32e·exp(4|T| + 4)
  (S2997: ‖E‖ ≤ (1+|t|)(16+12|t|)·2e^{π|t|/2} ≤ 32e^{4|t|}; |G_T| = e^{σ²−(t−T)²}; 4|t| − (t−T)² ≤ 4|T| + 4).
* `zetaE_gauss_bddAbove` — BddAbove ((‖·‖ ∘ F_T) '' verticalClosedStrip 0 1)  (compact part |s| ≤ 2 by continuity).
This discharges S2984's `hB`; with S2993 (edges, all t) + S2986/S2988 (Gaussian localisation) the convexity bound
|ζ(σ+iT)| ≪ (1+|T|)^{(1−σ)/2} log(2+|T|) is assembled next.
-/

open Complex Set Metric

namespace ZetaThreeLinesHB

theorem norm_F_far (T : ℝ) {s : ℂ} (h0 : 0 ≤ s.re) (h1 : s.re ≤ 1) (ht : 1 ≤ |s.im|) :
    ‖ZetaEntire.zetaE s * Complex.exp ((s - (T : ℂ) * I) ^ 2)‖
      ≤ 32 * Real.exp 1 * Real.exp (4 * |T| + 4) := by
  rw [norm_mul, ZetaThreeLines.norm_gauss]
  set u := |s.im| with hu
  have hE := ZetaEStrip.norm_zetaE_strip_le h0 h1 ht
  rw [← hu] at hE
  have hu0 : 0 ≤ u := abs_nonneg _
  have hpoly : (1 + u) * (16 + 12 * u) ≤ 16 * Real.exp (2 * u) := by
    have he : 1 + u ≤ Real.exp u := by linarith [Real.add_one_le_exp u]
    have h2 : (1 + u) ^ 2 ≤ Real.exp u ^ 2 := pow_le_pow_left₀ (by linarith) he 2
    rw [← Real.exp_nat_mul] at h2; push_cast at h2
    nlinarith
  have hpi : Real.pi * u / 2 ≤ 2 * u := by nlinarith [Real.pi_lt_four]
  have hE2 : ‖ZetaEntire.zetaE s‖ ≤ 32 * Real.exp (4 * u) := by
    refine le_trans hE ?_
    have h3 : Real.exp (Real.pi * u / 2) ≤ Real.exp (2 * u) := Real.exp_le_exp.mpr hpi
    calc (1 + u) * ((16 + 12 * u) * (2 * Real.exp (Real.pi * u / 2)))
        = 2 * ((1 + u) * (16 + 12 * u)) * Real.exp (Real.pi * u / 2) := by ring
      _ ≤ 2 * (16 * Real.exp (2 * u)) * Real.exp (2 * u) := by gcongr
      _ = 32 * Real.exp (4 * u) := by
          have e4 : Real.exp (4 * u) = Real.exp (2 * u) * Real.exp (2 * u) := by
            rw [← Real.exp_add]; ring_nf
          rw [e4]; ring
  have hσ : s.re ^ 2 ≤ 1 := by nlinarith
  have hg : Real.exp (s.re ^ 2 - (s.im - T) ^ 2) ≤ Real.exp 1 * Real.exp (-(s.im - T) ^ 2) := by
    rw [← Real.exp_add]; apply Real.exp_le_exp.mpr; linarith
  have hkey : Real.exp (4 * u) * Real.exp (-(s.im - T) ^ 2) ≤ Real.exp (4 * |T| + 4) := by
    rw [← Real.exp_add]; apply Real.exp_le_exp.mpr
    have htri : u ≤ |T| + |s.im - T| := by
      rw [hu]; have := abs_add_le T (s.im - T); simpa using this
    have hsq : (s.im - T) ^ 2 = |s.im - T| ^ 2 := (sq_abs _).symm
    nlinarith [sq_nonneg (|s.im - T| - 2)]
  calc ‖ZetaEntire.zetaE s‖ * Real.exp (s.re ^ 2 - (s.im - T) ^ 2)
      ≤ (32 * Real.exp (4 * u)) * (Real.exp 1 * Real.exp (-(s.im - T) ^ 2)) :=
        mul_le_mul hE2 hg (by positivity) (by positivity)
    _ = 32 * Real.exp 1 * (Real.exp (4 * u) * Real.exp (-(s.im - T) ^ 2)) := by ring
    _ ≤ 32 * Real.exp 1 * Real.exp (4 * |T| + 4) := by gcongr

theorem zetaE_gauss_bddAbove (T : ℝ) :
    BddAbove ((norm ∘ fun s => ZetaEntire.zetaE s * Complex.exp ((s - (T : ℂ) * I) ^ 2)) ''
      Complex.HadamardThreeLines.verticalClosedStrip 0 1) := by
  set F : ℂ → ℂ := fun s => ZetaEntire.zetaE s * Complex.exp ((s - (T : ℂ) * I) ^ 2) with hF
  have hc : Continuous F :=
    (ZetaEntire.zetaE_differentiable.mul
      ((((differentiable_id.sub_const _)).pow 2).cexp)).continuous
  obtain ⟨M, hM⟩ := (isCompact_closedBall (0 : ℂ) 2).exists_bound_of_continuousOn hc.continuousOn
  refine ⟨max M (32 * Real.exp 1 * Real.exp (4 * |T| + 4)), ?_⟩
  rintro _ ⟨s, hs, rfl⟩
  simp only [Complex.HadamardThreeLines.verticalClosedStrip, mem_preimage, mem_Icc] at hs
  simp only [Function.comp]
  rcases le_or_gt 1 |s.im| with ht | ht
  · exact le_trans (norm_F_far T hs.1 hs.2 ht) (le_max_right _ _)
  · have hball : s ∈ closedBall (0 : ℂ) 2 := by
      rw [mem_closedBall, dist_zero_right]
      have := Complex.norm_le_abs_re_add_abs_im s
      rw [abs_of_nonneg hs.1] at this
      linarith [hs.2]
    exact le_trans (hM s hball) (le_max_left _ _)

end ZetaThreeLinesHB

#print axioms ZetaThreeLinesHB.norm_F_far
#print axioms ZetaThreeLinesHB.zetaE_gauss_bddAbove
