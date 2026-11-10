import Mathlib

/-!
# S3034 — Schwarz-integral formula for the derivative (lane 1, Poisson–Landau program, brick P1)

For H analytic on ball 0 R' and 0 < R < R':
  H'(0) = circleAverage (ζ ↦ (2ζ/ζ²)·Re H(ζ)) over the circle of radius R,
i.e. H'(0) = (1/(πR))∫₀^{2π} Re H(Re^{iθ}) e^{−iθ} dθ — the derivative is determined by the REAL PART on the circle.
Proof: S(w) := circleAverage(HR(0,w,·)·Re H) (Herglotz–Riesz integral, Mathlib) is analytic with Re S = Re H on the disk
(harmonic Poisson formula, Re H harmonic); F = S − H has Re F ≡ 0 near 0, so F'(0) is both real (t ↦ F(t)) and imaginary
(t ↦ F(it)) ⇒ F'(0) = 0; S'(0) from `hasDerivAt_circleAverage_herglotzRieszKernel_smul`.
Use (Report): with H = log G, G = ζ/B zero-free (S3026), Re H = log|ζ| on the circle — so −Re ζ'/ζ is controlled by the
CIRCLE AVERAGE of log|ζ|·cos θ, not its maximum: left half ≤ convexity (1−Re z)/2·log t, right half ≥ −log(3+1/(σ−1)) ⇒ κ → 1/4.
-/

open Complex Metric Real

namespace SchwarzDeriv

theorem deriv_eq_zero_of_re_eq_zero {F : ℂ → ℂ} {F' : ℂ} {r : ℝ} (hr : 0 < r) (hF : HasDerivAt F F' 0)
    (hre : ∀ z ∈ ball (0 : ℂ) r, (F z).re = 0) : F' = 0 := by
  -- real direction
  have h1 : HasDerivAt (fun x : ℝ => (F x).re) F'.re 0 := by
    have : HasDerivAt F F' ((0 : ℝ) : ℂ) := by simpa using hF
    exact this.real_of_complex
  have hev1 : (fun x : ℝ => (F x).re) =ᶠ[nhds 0] fun _ => (0 : ℝ) := by
    have : Metric.ball (0 : ℝ) r ∈ nhds (0 : ℝ) := ball_mem_nhds 0 hr
    filter_upwards [this] with x hx
    apply hre
    simpa [mem_ball, dist_zero_right] using hx
  have hre0 : F'.re = 0 := (h1.congr_of_eventuallyEq hev1.symm).unique (hasDerivAt_const _ _) |>.symm ▸ rfl
  -- imaginary direction
  have hI : HasDerivAt (fun z : ℂ => F (I * z)) (F' * I) 0 := by
    have hlin : HasDerivAt (fun z : ℂ => I * z) I 0 := by simpa using (hasDerivAt_id (0 : ℂ)).const_mul I
    have hF0 : HasDerivAt F F' (I * 0) := by rw [mul_zero]; exact hF
    exact hF0.comp (0 : ℂ) hlin
  have h2 : HasDerivAt (fun x : ℝ => (F (I * x)).re) (F' * I).re 0 := by
    have : HasDerivAt (fun z : ℂ => F (I * z)) (F' * I) ((0 : ℝ) : ℂ) := by simpa using hI
    exact this.real_of_complex
  have hev2 : (fun x : ℝ => (F (I * x)).re) =ᶠ[nhds 0] fun _ => (0 : ℝ) := by
    have : Metric.ball (0 : ℝ) r ∈ nhds (0 : ℝ) := ball_mem_nhds 0 hr
    filter_upwards [this] with x hx
    apply hre
    simpa [mem_ball, dist_zero_right, norm_mul] using hx
  have him0 : (F' * I).re = 0 := (h2.congr_of_eventuallyEq hev2.symm).unique (hasDerivAt_const _ _) |>.symm ▸ rfl
  simp at him0
  exact Complex.ext hre0 him0

theorem deriv_eq_schwarz {H : ℂ → ℂ} {R R' : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hH : AnalyticOnNhd ℂ H (ball 0 R')) :
    deriv H 0 = circleAverage (fun ζ => (2 * ζ / (ζ - 0) ^ 2) • (((H ζ).re : ℝ) : ℂ)) 0 R := by
  set u : ℂ → ℝ := fun ζ => (H ζ).re with hu
  have hharm : InnerProductSpace.HarmonicOnNhd u (Metric.closedBall 0 R) := fun z hz =>
    (hH z (closedBall_subset_ball hRR' hz)).harmonicAt_re
  have hcont : ContinuousOn (fun ζ => ((u ζ : ℝ) : ℂ)) (sphere (0 : ℂ) R) := by
    have : ContinuousOn H (sphere (0 : ℂ) R) :=
      (hH.continuousOn).mono (sphere_subset_closedBall.trans (closedBall_subset_ball hRR'))
    exact (Complex.continuous_ofReal.comp_continuousOn (Complex.continuous_re.comp_continuousOn this))
  have hci : CircleIntegrable (fun ζ => ((u ζ : ℝ) : ℂ)) 0 R := hcont.circleIntegrable hR.le
  have hciu : CircleIntegrable u 0 R := by
    have : ContinuousOn u (sphere (0 : ℂ) R) :=
      Complex.continuous_re.comp_continuousOn
        ((hH.continuousOn).mono (sphere_subset_closedBall.trans (closedBall_subset_ball hRR')))
    exact this.circleIntegrable hR.le
  set S : ℂ → ℂ := fun w => circleAverage (fun ζ => herglotzRieszKernel 0 w ζ • ((u ζ : ℝ) : ℂ)) 0 R with hS
  have hSd : HasDerivAt S (circleAverage (fun ζ => (2 * ζ / (ζ - 0) ^ 2) • ((u ζ : ℝ) : ℂ)) 0 R) 0 :=
    hasDerivAt_circleAverage_herglotzRieszKernel_smul hci (mem_ball_self hR)
  have hSre : ∀ w ∈ ball (0 : ℂ) R, (S w).re = u w := by
    intro w hw
    rw [hS, re_circleAverage_herglotzRieszKernel_smul hciu hw]
    exact hharm.circleAverage_re_herglotzRieszKernel_smul hw
  have hHd : HasDerivAt H (deriv H 0) 0 :=
    (hH 0 (mem_ball_self (hR.trans hRR'))).differentiableAt.hasDerivAt
  have hF := hSd.sub hHd
  have hz := deriv_eq_zero_of_re_eq_zero hR hF (fun z hz => by
    simp only [Pi.sub_apply, Complex.sub_re, hSre z hz, hu, sub_self])
  rw [sub_eq_zero] at hz
  exact hz.symm

end SchwarzDeriv

#print axioms SchwarzDeriv.deriv_eq_zero_of_re_eq_zero
#print axioms SchwarzDeriv.deriv_eq_schwarz
