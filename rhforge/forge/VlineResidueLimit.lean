import VlineResidueShift

set_option maxHeartbeats 1000000

open Complex intervalIntegral Filter Topology

/-!
# §1661 ★★ — THE RESIDUE-DEPOSIT LIMIT: `∫_{Re=σ₀} c·(s-p)⁻¹ - ∫_{Re=σ₁} c·(s-p)⁻¹ → 2πc` as `T→∞`

The `T → ∞` limit of §1660: the horizontal connectors decay like `‖c‖/(T - |Im p|)` (uniform in the
integration variable, since the denominator's imaginary part dominates), so the truncated
vertical-line difference across the pole converges to the clean residue deposit `2πc`. With
`c = x²/2`, `p = 1` this is the Perron pole-term contour shift completed: moving the Perron line
from `σ₀ > 1` into the zero-free region `σ₁ < 1` deposits the PNT main term.
-/

theorem vline_diff_residue_tendsto (c p : ℂ) {σ₁ σ₀ : ℝ}
    (h1 : σ₁ < p.re) (h2 : p.re < σ₀) :
    Tendsto (fun T : ℝ =>
        (∫ y in (-T)..T, c * ((σ₀ : ℂ) + y * Complex.I - p)⁻¹) -
          (∫ y in (-T)..T, c * ((σ₁ : ℂ) + y * Complex.I - p)⁻¹))
      atTop (𝓝 (2 * (Real.pi : ℂ) * c)) := by
  -- the two horizontal connectors and their decay
  set bot : ℝ → ℂ := fun T => ∫ t in σ₁..σ₀, c * ((t : ℂ) + ((-T : ℝ) : ℂ) * Complex.I - p)⁻¹
    with hbotdef
  set top : ℝ → ℂ := fun T => ∫ t in σ₁..σ₀, c * ((t : ℂ) + ((T : ℝ) : ℂ) * Complex.I - p)⁻¹
    with htopdef
  have hg : Tendsto (fun T : ℝ => ‖c‖ / (T - |p.im|) * |σ₀ - σ₁|) atTop (𝓝 0) := by
    have hden : Tendsto (fun T : ℝ => T - |p.im|) atTop atTop := by
      simpa [sub_eq_add_neg] using tendsto_atTop_add_const_right atTop (-|p.im|) tendsto_id
    have hd := Tendsto.div_atTop (tendsto_const_nhds (x := ‖c‖)) hden
    simpa using hd.mul_const |σ₀ - σ₁|
  have htop0 : Tendsto top atTop (𝓝 0) := by
    apply squeeze_zero_norm'
    · filter_upwards [eventually_gt_atTop |p.im|] with T hT
      refine intervalIntegral.norm_integral_le_of_norm_le_const
        (C := ‖c‖ / (T - |p.im|)) fun t _ => ?_
      have him : ((t : ℂ) + ((T : ℝ) : ℂ) * Complex.I - p).im = T - p.im := by simp
      have hz : T - |p.im| ≤ ‖(t : ℂ) + ((T : ℝ) : ℂ) * Complex.I - p‖ := by
        have habs := Complex.abs_im_le_norm ((t : ℂ) + ((T : ℝ) : ℂ) * Complex.I - p)
        rw [him] at habs
        have h5 := le_abs_self (T - p.im)
        have h6 := le_abs_self p.im
        linarith
      rw [norm_mul, norm_inv, div_eq_mul_inv]
      apply mul_le_mul_of_nonneg_left _ (norm_nonneg c)
      rw [inv_eq_one_div, inv_eq_one_div]
      exact one_div_le_one_div_of_le (sub_pos.mpr hT) hz
    · exact hg
  have hbot0 : Tendsto bot atTop (𝓝 0) := by
    apply squeeze_zero_norm'
    · filter_upwards [eventually_gt_atTop |p.im|] with T hT
      refine intervalIntegral.norm_integral_le_of_norm_le_const
        (C := ‖c‖ / (T - |p.im|)) fun t _ => ?_
      have him : ((t : ℂ) + ((-T : ℝ) : ℂ) * Complex.I - p).im = -T - p.im := by simp
      have hz : T - |p.im| ≤ ‖(t : ℂ) + ((-T : ℝ) : ℂ) * Complex.I - p‖ := by
        have habs := Complex.abs_im_le_norm ((t : ℂ) + ((-T : ℝ) : ℂ) * Complex.I - p)
        rw [him] at habs
        have h6 := neg_abs_le p.im
        have h7 : -T - p.im < 0 := by linarith [neg_abs_le p.im, le_abs_self p.im]
        rw [abs_of_neg h7] at habs
        linarith
      rw [norm_mul, norm_inv, div_eq_mul_inv]
      apply mul_le_mul_of_nonneg_left _ (norm_nonneg c)
      rw [inv_eq_one_div, inv_eq_one_div]
      exact one_div_le_one_div_of_le (sub_pos.mpr hT) hz
    · exact hg
  -- eventual identity from §1660
  have key : ∀ᶠ T : ℝ in atTop,
      2 * (Real.pi : ℂ) * c + Complex.I * (bot T - top T) =
        (∫ y in (-T)..T, c * ((σ₀ : ℂ) + y * Complex.I - p)⁻¹) -
          (∫ y in (-T)..T, c * ((σ₁ : ℂ) + y * Complex.I - p)⁻¹) := by
    filter_upwards [eventually_gt_atTop |p.im|] with T hT
    have ha : -T < p.im := by linarith [neg_abs_le p.im]
    have hb : p.im < T := by linarith [le_abs_self p.im]
    exact (vline_diff_residue c p h1 h2 ha hb).symm
  -- assemble
  have hlim : Tendsto (fun T : ℝ => 2 * (Real.pi : ℂ) * c + Complex.I * (bot T - top T))
      atTop (𝓝 (2 * (Real.pi : ℂ) * c)) := by
    have h0 : Tendsto (fun T : ℝ => Complex.I * (bot T - top T)) atTop (𝓝 0) := by
      have := (hbot0.sub htop0).const_mul Complex.I
      simpa using this
    simpa using tendsto_const_nhds.add h0
  exact hlim.congr' key

#print axioms vline_diff_residue_tendsto
