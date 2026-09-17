import Mathlib.Analysis.Complex.JensenFormula

/-!
# §1778 — Jensen no-interior-zeros criterion (the analytic bridge for the CS engine)

§1777 (`SelfInversiveOnCircle`) reduced the Connes–van Suijlekom / Carathéodory-Fejér on-circle
conclusion to ONE owed input: the (self-inversive) kernel polynomial has NO zeros strictly inside the
open unit disk.  This brick supplies the analytic criterion that discharges it, built directly on
Mathlib's Jensen formula (`AnalyticOnNhd.circleAverage_log_norm`):

    circleAverage (log‖f‖) 0 1  =  Σᶠ u, divisor(u) · log(1 · ‖0 - u‖⁻¹)  +  log‖f 0‖.

Every term of the divisor sum is `≥ 0` (`divisor ≥ 0` for an analytic `f`, and `log(1/‖u‖) ≥ 0` for
`‖u‖ ≤ 1`), so the boundary log-modulus average always DOMINATES `log‖f 0‖` by the "Jensen excess" =
the interior zero-counting sum.  Hence if that average has NO excess (`≤ log‖f 0‖`), the divisor must
vanish at every point of the open disk: **no interior zeros**.  This is exactly the Szegő /
positive-measure condition that the PSD-singular Toeplitz structure delivers (the Carathéodory
content), now reduced to a boundary log-modulus equality — the analytic half feeding §1777.
EXIT0 / std-axioms.
-/

open Complex Real MeasureTheory Metric

/-- **Jensen no-interior-zeros criterion.**  If `f` is analytic on the closed unit disk with
`f 0 ≠ 0`, and the circle-average of `log‖f‖` over the unit circle does not exceed `log‖f 0‖`
(no Jensen excess), then the divisor of `f` vanishes at every point of the open unit disk. -/
theorem jensen_no_interior_zeros {f : ℂ → ℂ}
    (h₁f : AnalyticOnNhd ℂ f (closedBall 0 1))
    (h₂f : f 0 ≠ 0)
    (havg : circleAverage (fun z => Real.log ‖f z‖) 0 1 ≤ Real.log ‖f 0‖)
    (u : ℂ) (hu : ‖u‖ < 1) :
    MeromorphicOn.divisor f (closedBall 0 1) u = 0 := by
  set D := MeromorphicOn.divisor f (closedBall 0 1) with hD
  -- Jensen's formula (c = 0, R = 1); normalize `|1|` to `1`
  have h₁f' : AnalyticOnNhd ℂ f (closedBall 0 |(1 : ℝ)|) := by rwa [abs_one]
  have hJ := h₁f'.circleAverage_log_norm (c := 0) (R := 1) (by norm_num) h₂f
  rw [abs_one] at hJ
  rw [← hD] at hJ
  -- rewrite the no-excess hypothesis through Jensen
  rw [hJ] at havg
  -- the divisor-sum summand
  set g : ℂ → ℝ := fun v => (D v : ℝ) * Real.log (1 * ‖(0 : ℂ) - v‖⁻¹) with hg
  have hDnonneg : 0 ≤ D := MeromorphicOn.AnalyticOnNhd.divisor_nonneg h₁f
  -- each summand is nonneg
  have hgnonneg : ∀ v, 0 ≤ g v := by
    intro v
    simp only [hg]
    by_cases hDv0 : D v = 0
    · simp [hDv0]
    · have hvmem : v ∈ closedBall (0 : ℂ) 1 :=
        D.supportWithinDomain hDv0
      have hnorm : ‖v‖ ≤ 1 := mem_closedBall_zero_iff.mp hvmem
      have hDv : 0 ≤ (D v : ℝ) := by exact_mod_cast hDnonneg v
      have hlog : 0 ≤ Real.log (1 * ‖(0 : ℂ) - v‖⁻¹) := by
        rw [one_mul, zero_sub, norm_neg, Real.log_inv]
        exact neg_nonneg.mpr (Real.log_nonpos (norm_nonneg v) hnorm)
      positivity
  -- finite support of g (⊆ support of D, finite on the compact ball)
  have hgfin : (Function.support g).Finite := by
    apply (D.finiteSupport (isCompact_closedBall 0 1)).subset
    intro v hv
    simp only [hg, Function.mem_support, ne_eq, mul_eq_zero, not_or] at hv
    exact Function.mem_support.mpr (Int.cast_ne_zero.mp hv.1)
  -- from Jensen + no-excess: the divisor sum is ≤ 0, hence = 0
  have hsum_nonneg : 0 ≤ ∑ᶠ v, g v := finsum_nonneg hgnonneg
  have hsum_le : (∑ᶠ v, g v) ≤ 0 := by linarith [havg]
  have hsum0 : (∑ᶠ v, g v) = 0 := le_antisymm hsum_le hsum_nonneg
  -- single term ≤ total = 0
  have hgu_le : g u ≤ 0 := hsum0 ▸ single_le_finsum u hgfin hgnonneg
  have hDu_nonneg : 0 ≤ (D u : ℝ) := by exact_mod_cast hDnonneg u
  -- split on whether u = 0
  by_cases hu0 : u = 0
  · subst hu0
    rw [hD, MeromorphicOn.AnalyticOnNhd.divisor_apply h₁f (by simp : (0 : ℂ) ∈ closedBall (0 : ℂ) 1)]
    simp [(h₁f 0 (by simp)).analyticOrderAt_eq_zero.mpr h₂f]
  · have hupos : 0 < ‖u‖ := norm_pos_iff.mpr hu0
    have hlogu_pos : 0 < Real.log (1 * ‖(0 : ℂ) - u‖⁻¹) := by
      rw [one_mul, zero_sub, norm_neg, Real.log_inv]
      exact neg_pos.mpr (Real.log_neg hupos hu)
    have hDu_le : (D u : ℝ) ≤ 0 := by
      by_contra h
      push_neg at h
      have : 0 < g u := by simp only [hg]; positivity
      linarith
    have hcast : (D u : ℝ) = 0 := le_antisymm hDu_le hDu_nonneg
    exact_mod_cast hcast
