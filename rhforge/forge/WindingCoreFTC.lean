import Mathlib.Analysis.SpecialFunctions.Trigonometric.ArctanDeriv
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Tactic

set_option maxHeartbeats 1000000

open Complex intervalIntegral

/-!
# §1654 — WINDING CORE: closed form for `∫ (t + k·I)⁻¹ dt` (k ≠ 0), branch-cut-free

The base integral of the rectangle-winding arc. Mathlib has NO rectangle-winding lemma (only circle
winding + rectangle Cauchy-Goursat VANISHING); this arc builds it. The key idea: avoid `Complex.log`
branch cuts entirely by integrating the REAL decomposition `(t+ki)⁻¹ = (t - ki)/(t²+k²)` — real part
antiderivative `½log(t²+k²)`, imaginary part antiderivative `-arctan(t/k)` — both globally smooth
for `k ≠ 0`. All four rectangle edges reduce to this one closed form.
-/

/-- The explicit antiderivative of `t ↦ ((t:ℂ) + k·I)⁻¹`. -/
noncomputable def windingPrim (k t : ℝ) : ℂ :=
  (↑(Real.log (t ^ 2 + k ^ 2) / 2) : ℂ) - Complex.I * ↑(Real.arctan (t / k))

theorem windingLine_ne_zero {k : ℝ} (hk : k ≠ 0) (t : ℝ) : (t : ℂ) + k * Complex.I ≠ 0 := by
  intro h
  have him := congrArg Complex.im h
  simp at him
  exact hk him

theorem hasDerivAt_windingPrim {k : ℝ} (hk : k ≠ 0) (t : ℝ) :
    HasDerivAt (windingPrim k) (((t : ℂ) + k * Complex.I)⁻¹) t := by
  have hpos : (0 : ℝ) < t ^ 2 + k ^ 2 := by positivity
  -- real part: d/dt [log(t²+k²)/2] = t/(t²+k²)
  have hpoly : HasDerivAt (fun s : ℝ => s ^ 2 + k ^ 2) (2 * t) t := by
    simpa using (hasDerivAt_pow 2 t).add_const (k ^ 2)
  have hlog : HasDerivAt (fun s : ℝ => Real.log (s ^ 2 + k ^ 2) / 2) (t / (t ^ 2 + k ^ 2)) t := by
    have h1 := (Real.hasDerivAt_log hpos.ne').comp t hpoly
    have h2 := h1.div_const 2
    have hval : (t ^ 2 + k ^ 2)⁻¹ * (2 * t) / 2 = t / (t ^ 2 + k ^ 2) := by
      ring
    rwa [hval] at h2
  -- imaginary part: d/dt [arctan(t/k)] = k/(t²+k²)
  have harc : HasDerivAt (fun s : ℝ => Real.arctan (s / k)) (k / (t ^ 2 + k ^ 2)) t := by
    have hin : HasDerivAt (fun s : ℝ => s / k) (1 / k) t := by
      simpa using (hasDerivAt_id t).div_const k
    have h1 := (Real.hasDerivAt_arctan (t / k)).comp t hin
    have hval : 1 / (1 + (t / k) ^ 2) * (1 / k) = k / (t ^ 2 + k ^ 2) := by
      have hk2 : k ^ 2 ≠ 0 := pow_ne_zero 2 hk
      field_simp
      ring
    rwa [hval] at h1
  -- lift to ℂ and combine
  have hlogC := hlog.ofReal_comp
  have harcC := harc.ofReal_comp
  have hG := hlogC.sub (harcC.const_mul Complex.I)
  have hden : ((t : ℂ) ^ 2 + (k : ℂ) ^ 2) ≠ 0 := by
    have : ((t ^ 2 + k ^ 2 : ℝ) : ℂ) ≠ 0 := by exact_mod_cast hpos.ne'
    push_cast at this
    exact this
  have hprod : ((t : ℂ) + k * Complex.I) * (((t : ℂ) - k * Complex.I) / ((t : ℂ) ^ 2 + (k : ℂ) ^ 2))
      = 1 := by
    rw [mul_div_assoc', div_eq_one_iff_eq hden]
    linear_combination (-(k : ℂ) ^ 2) * Complex.I_sq
  have hval : (↑(t / (t ^ 2 + k ^ 2)) : ℂ) - Complex.I * ↑(k / (t ^ 2 + k ^ 2))
      = ((t : ℂ) + k * Complex.I)⁻¹ := by
    rw [inv_eq_of_mul_eq_one_right hprod]
    push_cast
    ring
  rw [← hval]
  exact hG

/-- §1654 — THE CORE CLOSED FORM: for `k ≠ 0`,
`∫_{t₁}^{t₂} (t + k·I)⁻¹ dt = ½log((t₂²+k²)/(t₁²+k²)) - i(arctan(t₂/k) - arctan(t₁/k))`
(stated via `windingPrim` differences; completely branch-cut-free). -/
theorem integral_winding_line {k : ℝ} (hk : k ≠ 0) (t₁ t₂ : ℝ) :
    ∫ t in t₁..t₂, ((t : ℂ) + k * Complex.I)⁻¹ = windingPrim k t₂ - windingPrim k t₁ := by
  have hcont : Continuous fun t : ℝ => ((t : ℂ) + k * Complex.I)⁻¹ := by
    apply Continuous.inv₀
    · exact Complex.continuous_ofReal.add continuous_const
    · exact windingLine_ne_zero hk
  exact integral_eq_sub_of_hasDerivAt (fun t _ => hasDerivAt_windingPrim hk t)
    (hcont.intervalIntegrable t₁ t₂)

#print axioms integral_winding_line
