import CosHalfIntegrals

/-!
# S3041 — the circle-average bound: −avg Re((2ζ/ζ²)·u) ≤ L/4 + 2(A+D)/(πr)  (lane 1, Poisson–Landau brick P3)

For u continuous on the circle |ζ| = r with
  u ζ ≤ A + (L/2)·max(−Re ζ, 0)      (left-half growth: convexity, slope L/2 in the depth −Re ζ)
  Re ζ ≥ 0 ⇒ −D ≤ u ζ                (right-half floor: |ζ(s)| ≥ (3 + 1/(σ−1))⁻¹ for Re s > 1)
we get −avg_{|ζ|=r} Re((2ζ/ζ²)·u(ζ)) ≤ L/4 + 2(A + D)/(π r).
Pointwise Re((2ζ/ζ²)·u) = (2cos θ/r)·u; split by the sign of cos θ; integrate with S3036
(∫max(−cos,0)² = π/2, ∫max(±cos,0) = 2).  With L = log t this is the κ = 1/4 growth term.
-/

open Complex Real intervalIntegral MeasureTheory Metric

namespace PoissonAvg

lemma re_kernel_circleMap {r : ℝ} (hr : 0 < r) (θ : ℝ) (x : ℝ) :
    ((2 * circleMap 0 r θ / (circleMap 0 r θ - 0) ^ 2) * ((x : ℝ) : ℂ)).re = 2 * Real.cos θ / r * x := by
  set ζ := circleMap 0 r θ with hζ
  have hne : ζ ≠ 0 := by simp [hζ, circleMap, hr.ne', Complex.exp_ne_zero]
  have hn : Complex.normSq ζ = r ^ 2 := by
    rw [Complex.normSq_eq_norm_sq, hζ, norm_circleMap_zero, abs_of_pos hr]
  have hre : ζ.re = r * Real.cos θ := by
    simp [hζ, circleMap, Complex.mul_re, ← Complex.ofReal_cos, ← Complex.ofReal_sin, Complex.exp_mul_I]
  rw [sub_zero, show 2 * ζ / ζ ^ 2 = 2 / ζ by field_simp, Complex.re_mul_ofReal, Complex.div_re]
  simp only [hn, hre]
  norm_num
  field_simp
  try simp

theorem neg_avg_le {u : ℂ → ℝ} {r A D L : ℝ} (hr : 0 < r) (hA : 0 ≤ A) (hD : 0 ≤ D) (hL : 0 ≤ L)
    (hcont : ContinuousOn u (sphere (0 : ℂ) r))
    (hup : ∀ ζ ∈ sphere (0 : ℂ) r, u ζ ≤ A + L / 2 * max (-ζ.re) 0)
    (hlo : ∀ ζ ∈ sphere (0 : ℂ) r, 0 ≤ ζ.re → -D ≤ u ζ) :
    -Real.circleAverage (fun ζ => ((2 * ζ / (ζ - 0) ^ 2) * ((u ζ : ℝ) : ℂ)).re) 0 r
      ≤ L / 4 + 2 * (A + D) / (π * r) := by
  rw [Real.circleAverage_def]
  set U : ℝ → ℝ := fun θ => u (circleMap 0 r θ) with hU
  have hmem : ∀ θ, circleMap 0 r θ ∈ sphere (0 : ℂ) r := fun θ => by
    simpa [abs_of_pos hr] using circleMap_mem_sphere' (0 : ℂ) r θ
  have hre : ∀ θ, (circleMap 0 r θ).re = r * Real.cos θ := fun θ => by
    simp [circleMap, Complex.mul_re, ← Complex.ofReal_cos, ← Complex.ofReal_sin, Complex.exp_mul_I]
  have hUc : Continuous U := hcont.comp_continuous (continuous_circleMap 0 r) hmem
  have hint : (fun θ => ((2 * circleMap 0 r θ / (circleMap 0 r θ - 0) ^ 2) * ((u (circleMap 0 r θ) : ℝ) : ℂ)).re)
      = fun θ => 2 * Real.cos θ / r * U θ := by
    funext θ; exact re_kernel_circleMap hr θ _
  simp only [Function.comp_def, smul_eq_mul]
  rw [hint]
  set Bd : ℝ → ℝ := fun θ => 2 / r * (A * max (-Real.cos θ) 0 + L * r / 2 * (max (-Real.cos θ) 0) ^ 2
      + D * max (Real.cos θ) 0) with hBd
  have hpt : ∀ θ, -(2 * Real.cos θ / r * U θ) ≤ Bd θ := by
    intro θ
    have h2r : 0 < 2 / r := by positivity
    have key : -(Real.cos θ * U θ) ≤ A * max (-Real.cos θ) 0 + L * r / 2 * (max (-Real.cos θ) 0) ^ 2
        + D * max (Real.cos θ) 0 := by
      rcases le_or_gt (Real.cos θ) 0 with hc | hc
      · have hu := hup _ (hmem θ)
        rw [hre θ] at hu
        rw [max_eq_left (by linarith : (0 : ℝ) ≤ -Real.cos θ), max_eq_right hc] at *
        have : max (-(r * Real.cos θ)) 0 = r * (-Real.cos θ) := by
          rw [max_eq_left (by nlinarith)]; ring
        rw [this] at hu
        nlinarith
      · have hu := hlo _ (hmem θ) (by rw [hre θ]; positivity)
        rw [max_eq_right (by linarith : -Real.cos θ ≤ 0), max_eq_left hc.le]
        nlinarith
    have : -(2 * Real.cos θ / r * U θ) = 2 / r * (-(Real.cos θ * U θ)) := by ring
    rw [this, hBd]
    exact mul_le_mul_of_nonneg_left key h2r.le
  have hcf : Continuous fun θ => 2 * Real.cos θ / r * U θ := by fun_prop
  have hBc : Continuous Bd := by rw [hBd]; fun_prop
  have hmono : ∫ θ in (0 : ℝ)..(2 * π), -(2 * Real.cos θ / r * U θ) ≤ ∫ θ in (0 : ℝ)..(2 * π), Bd θ :=
    integral_mono_on (by positivity) (hcf.neg.intervalIntegrable _ _) (hBc.intervalIntegrable _ _)
      (fun θ _ => hpt θ)
  have hBint : ∫ θ in (0 : ℝ)..(2 * π), Bd θ = 2 / r * (A * 2 + L * r / 2 * (π / 2) + D * 2) := by
    have e : ∫ θ in (0 : ℝ)..(2 * π), Bd θ = 2 / r * ∫ θ in (0 : ℝ)..(2 * π),
        (A * max (-Real.cos θ) 0 + L * r / 2 * (max (-Real.cos θ) 0) ^ 2 + D * max (Real.cos θ) 0) :=
      intervalIntegral.integral_const_mul _ _
    rw [e]
    congr 1
    have c1 : Continuous fun θ : ℝ => A * max (-Real.cos θ) 0 := by fun_prop
    have c2 : Continuous fun θ : ℝ => L * r / 2 * (max (-Real.cos θ) 0) ^ 2 := by fun_prop
    have c3 : Continuous fun θ : ℝ => D * max (Real.cos θ) 0 := by fun_prop
    rw [intervalIntegral.integral_add ((by fun_prop : Continuous fun θ : ℝ => A * max (-Real.cos θ) 0 + L * r / 2 * (max (-Real.cos θ) 0) ^ 2).intervalIntegrable _ _) (c3.intervalIntegrable _ _),
      intervalIntegral.integral_add (c1.intervalIntegrable _ _) (c2.intervalIntegrable _ _),
      intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
      CosHalf.integral_negcos_pos, CosHalf.integral_negcos_pos_sq, CosHalf.integral_cos_pos]
  rw [intervalIntegral.integral_neg] at hmono
  have hpi : 0 < π := pi_pos
  have h2pi : 0 < (2 * π)⁻¹ := by positivity
  have := mul_le_mul_of_nonneg_left hmono h2pi.le
  rw [hBint] at this
  have eq1 : -((2 * π)⁻¹ * ∫ θ in (0 : ℝ)..(2 * π), 2 * Real.cos θ / r * U θ)
      = (2 * π)⁻¹ * -∫ θ in (0 : ℝ)..(2 * π), 2 * Real.cos θ / r * U θ := by ring
  rw [eq1]
  have eq2 : (2 * π)⁻¹ * (2 / r * (A * 2 + L * r / 2 * (π / 2) + D * 2)) = L / 4 + 2 * (A + D) / (π * r) := by
    field_simp; ring
  linarith

end PoissonAvg

#print axioms PoissonAvg.re_kernel_circleMap
#print axioms PoissonAvg.neg_avg_le
