import NormLogDerivBallCenter
import NormFactorLe
import NormProdHalfBall
import Mathlib.Analysis.Complex.AbsMax
import Mathlib.NumberTheory.LSeries.RiemannZeta

set_option maxHeartbeats 1000000

open Complex Metric

/-!
# §1679 — remainder-norm bound, QUARTER-BALL geometry (the ζ-compatible parameterization)

§1672 required `g` nonvanishing on the half-ball while the zero-count (§1673/§1674) locates `S` in
the quarter-ball — for `ζ` the annulus `(R/4, R/2)` may contain zeros, breaking §1672's hypothesis.
This variant is the CONSISTENT parameterization: `S` = zeros in `closedBall c (R/4)`, `g`
nonvanishing only on `ball c (R/4)`, branch/BC/Cauchy run inside the quarter-ball, conclusion on
`closedBall c (R/16)` with constant `32/R`. Sphere division now gains the 3/4-gap:
`M' = Mζ/(3R/4)^Σm`.
-/

theorem norm_logDeriv_g_le_of_quarterBall_factor {c : ℂ} {R Mζ : ℝ} {g P : ℂ → ℂ} {S : Finset ℂ}
    {m : ℂ → ℕ} (hR : 0 < R) (hPdef : P = fun s => ∏ ρ ∈ S, (s - ρ) ^ (m ρ))
    (hS : ∀ ρ ∈ S, ρ ∈ closedBall c (R / 4))
    (hfac : ∀ s ∈ sphere c R, riemannZeta s = P s * g s)
    (hζsphere : ∀ s ∈ sphere c R, ‖riemannZeta s‖ ≤ Mζ)
    (hgdiff : DiffContOnCl ℂ g (ball c R)) (hgnz : ∀ z ∈ ball c (R / 4), g z ≠ 0)
    (hgc_lt : ‖g c‖ < Mζ / (3 * R / 4) ^ (∑ ρ ∈ S, m ρ)) :
    ∀ s ∈ closedBall c (R / 16),
      ‖logDeriv g s‖
        ≤ 32 * (Real.log (Mζ / (3 * R / 4) ^ (∑ ρ ∈ S, m ρ)) - Real.log ‖g c‖) / R := by
  set M' : ℝ := Mζ / (3 * R / 4) ^ (∑ ρ ∈ S, m ρ) with hM'def
  have hR4 : 0 < R / 4 := by linarith
  have h34 : (0 : ℝ) < 3 * R / 4 := by linarith
  have hgc0 : g c ≠ 0 := hgnz c (mem_ball_self hR4)
  have hgcpos : 0 < ‖g c‖ := norm_pos_iff.mpr hgc0
  have hM'pos : 0 < M' := lt_trans hgcpos hgc_lt
  -- sphere bound for g via the 3/4-gap
  have hgsphere : ∀ s ∈ sphere c R, ‖g s‖ ≤ M' := by
    intro s hs
    have hPlow : (3 * R / 4) ^ (∑ ρ ∈ S, m ρ) ≤ ‖P s‖ := by
      rw [hPdef, norm_prod]
      simp_rw [norm_pow]
      rw [← Finset.prod_pow_eq_pow_sum]
      apply Finset.prod_le_prod
      · intro ρ _
        positivity
      · intro ρ hρ
        have hbase : 3 * R / 4 ≤ ‖s - ρ‖ := by
          have hρc : dist ρ c ≤ R / 4 := mem_closedBall.mp (hS ρ hρ)
          have hsc : dist s c = R := mem_sphere.mp hs
          have htri := dist_triangle s ρ c
          rw [← dist_eq_norm]
          linarith
        gcongr
    exact norm_factor_le (hfac s hs) (by positivity) hPlow (hζsphere s hs)
  -- max modulus to the quarter-ball
  have hgball : ∀ z ∈ ball c (R / 4), ‖g z‖ ≤ M' := by
    intro z hz
    have hfr : ∀ u ∈ frontier (ball c R), ‖g u‖ ≤ M' := by
      intro u hu
      rw [frontier_ball c hR.ne'] at hu
      exact hgsphere u hu
    exact Complex.norm_le_of_forall_mem_frontier_norm_le isBounded_ball hgdiff hfr
      (subset_closure (by
        rw [mem_ball] at hz ⊢
        linarith))
  -- norm-Landau core on the quarter-ball
  set M'' : ℝ := Real.log (M' / ‖g c‖) with hM''def
  have hM''pos : 0 < M'' := Real.log_pos ((one_lt_div hgcpos).mpr hgc_lt)
  have hgrow : ∀ z ∈ ball c (R / 4), ‖g z‖ ≤ Real.exp M'' * ‖g c‖ := by
    intro z hz
    rw [hM''def, Real.exp_log (by positivity), div_mul_cancel₀ _ hgcpos.ne']
    exact hgball z hz
  have hgdiff2 : DifferentiableOn ℂ g (ball c (R / 4)) := by
    apply hgdiff.differentiableOn.mono
    intro z hz
    rw [mem_ball] at hz ⊢
    linarith
  have hcore := norm_logDeriv_le_of_forall_norm_le_center hR4 hM''pos hgdiff2 hgnz hgrow
  intro s hs
  have hs' : s ∈ closedBall c (R / 4 / 4) := by
    rw [mem_closedBall] at hs ⊢
    linarith
  calc ‖logDeriv g s‖ ≤ 8 * M'' / (R / 4) := hcore s hs'
    _ = 32 * M'' / R := by field_simp; ring
    _ = 32 * (Real.log M' - Real.log ‖g c‖) / R := by
        rw [hM''def, Real.log_div hM'pos.ne' hgcpos.ne']

#print axioms norm_logDeriv_g_le_of_quarterBall_factor
