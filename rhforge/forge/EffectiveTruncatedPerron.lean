import PerronVlineShift
import LeftLineIntegralParam
import LeftLinePackageParam
import ConnectorIntegralParam
import ConnectorIntegralBottom
import LogDerivZetaSmallT

set_option maxHeartbeats 4000000

open Complex MeasureTheory Set Metric

open scoped Interval

/-!
# §1698 ★★★★★ — THE EFFECTIVE TRUNCATED PERRON FORMULA

THE assembled theorem of the riemann9 overnight: there are `c₁, K > 0` and a threshold `T₂` such
that for all `T ≥ T₂`, `x ≥ 1`, `σ₀ ∈ (1, 2]`:

  `‖∫_{|t|≤T} F(σ₀+it) dt − 2π·(x²/2)‖ ≤ K·((1+log(T+2))·x^{σ₁(T)+1} + log²(T+2)·x³/T²)`

with `σ₁(T) = 1 - c₁/log(T+2)` and `F = -ζ'/ζ(s)·x^{s+1}/(s(s+1))` the smoothed-Perron integrand.
The truncated Perron integral IS the smoothed main term `2π·x²/2` up to the classical error: the
first term is `x²·x^{-c₁/log(T+2)}·O(log T)` — choosing `log T ≍ √log x` makes it
`x²·e^{-c√log x}` — and the second vanishes as `T → ∞`. Composition: the §1666 shift identity on
the zero-free box (§1697b package supplies zero-freeness), the §1697a left-line error integral, and
the §1697c/§1698a connector integrals, all with ONE coherent constant `c₁ = min(c_conn, c_landau)`.
Every ingredient of the chain was machine-checked tonight, ending here.
-/

theorem effective_truncated_perron :
    ∃ c₁ K T₂ : ℝ, 0 < c₁ ∧ 0 < K ∧ 6 ≤ T₂ ∧ ∀ T : ℝ, T₂ ≤ T → ∀ x : ℝ, 1 ≤ x →
      ∀ σ₀ : ℝ, 1 < σ₀ → σ₀ ≤ 2 →
      ‖(∫ t in (-T)..T,
          -logDeriv riemannZeta ((σ₀ : ℝ) + (t : ℝ) * Complex.I) *
              (x : ℂ) ^ (((σ₀ : ℝ) + (t : ℝ) * Complex.I) + 1) /
            (((σ₀ : ℝ) + (t : ℝ) * Complex.I) * (((σ₀ : ℝ) + (t : ℝ) * Complex.I) + 1)))
        - 2 * (Real.pi : ℂ) * ((x : ℂ) ^ (2 : ℂ) / 2)‖
      ≤ K * ((1 + Real.log (T + 2)) * x ^ (1 - c₁ / Real.log (T + 2) + 1)
          + (Real.log (T + 2)) ^ 2 * x ^ (3 : ℝ) / T ^ 2) := by
  obtain ⟨cA, A, T₀, hcA, hA, hT₀, hconn⟩ := connector_pointwise_bound
  obtain ⟨cB, hcB, hdist⟩ := landau_ball_zero_distance
  obtain ⟨σ₂, C, hσ₂half, hσ₂lt, hC, hsmall⟩ := logDeriv_zeta_small_t_bound T₀ (by linarith)
  set c₁ : ℝ := min cA cB with hc₁def
  have hc₁ : 0 < c₁ := lt_min hcA hcB
  have hc₁A : c₁ ≤ cA := min_le_left _ _
  have hpkg := left_line_package_param hcA hA hT₀ hconn hcB hdist hσ₂half hσ₂lt hC hsmall
  set T₁ : ℝ := max T₀ (Real.exp (min cA cB / (1 - σ₂)) + 6) with hT₁def
  have hT₁T₀ : T₀ ≤ T₁ := le_max_left _ _
  have hT₁6 : (6 : ℝ) ≤ T₁ := le_trans hT₀ hT₁T₀
  -- the coherent pointwise connector bound for c₁
  have hptmin : ∀ T : ℝ, T₀ ≤ T → ∀ σ : ℝ, 1 - c₁ / Real.log (T + 2) ≤ σ → σ ≤ 2 →
      ‖logDeriv riemannZeta ((σ : ℝ) + (T : ℝ) * Complex.I)‖
        ≤ A * (Real.log (T + 2)) ^ 2 := by
    intro T hT σ hσ h2
    apply hconn T hT σ _ h2
    have hL : 0 < Real.log (T + 2) := Real.log_pos (by linarith)
    have : c₁ / Real.log (T + 2) ≤ cA / Real.log (T + 2) := by gcongr
    linarith
  have hlineb := left_line_integral_bound_param hc₁ hA hT₁6 hT₀ hT₁T₀ hC hpkg
  have hconT := connector_integral_bound_param hc₁ hA hT₀ hptmin
  have hconB := connector_integral_bound_bottom hc₁ hA hT₀ hptmin
  -- the threshold
  refine ⟨c₁, (192 * A + 4 * T₀ * C + 4 * T₀ / c₁) + 4 * A,
    max (max T₁ (Real.exp (2 * c₁))) T₀, hc₁, by positivity,
    le_trans hT₁6 (le_trans (le_max_left _ _) (le_max_left _ _)), ?_⟩
  intro T hT x hx σ₀ hσ₀1 hσ₀2
  have hTT₁e : max T₁ (Real.exp (2 * c₁)) ≤ T := le_trans (le_max_left _ _) hT
  have hTT₀ : T₀ ≤ T := le_trans (le_max_right _ _) hT
  have hTT₁ : T₁ ≤ T := le_trans (le_max_left _ _) hTT₁e
  have hT6 : (6 : ℝ) ≤ T := le_trans hT₁6 hTT₁
  have hL : 0 < Real.log (T + 2) := Real.log_pos (by linarith)
  have hL2c : 2 * c₁ ≤ Real.log (T + 2) := by
    rw [← Real.log_exp (2 * c₁)]
    apply Real.log_le_log (Real.exp_pos _)
    have h := le_trans (le_max_right T₁ (Real.exp (2 * c₁))) hTT₁e
    linarith
  set L : ℝ := Real.log (T + 2) with hLdef
  set σ₁ : ℝ := 1 - c₁ / L with hσ₁def
  have hσ₁half : 1 / 2 ≤ σ₁ := by
    rw [hσ₁def]
    have : c₁ / L ≤ 1 / 2 := by
      rw [div_le_iff₀ hL]
      linarith
    linarith
  have hσ₁lt : σ₁ < 1 := by
    rw [hσ₁def]
    have : 0 < c₁ / L := by positivity
    linarith
  obtain ⟨_, hzf⟩ := hpkg T hTT₁
  -- the shift identity (§1666) on the zero-free box
  set z : ℂ := ((σ₁ : ℝ) : ℂ) + ((-T : ℝ) : ℂ) * Complex.I with hzdef
  set w : ℂ := ((σ₀ : ℝ) : ℂ) + ((T : ℝ) : ℂ) * Complex.I with hwdef
  have hz_re : z.re = σ₁ := by rw [hzdef]; simp
  have hz_im : z.im = -T := by rw [hzdef]; simp
  have hw_re : w.re = σ₀ := by rw [hwdef]; simp
  have hw_im : w.im = T := by rw [hwdef]; simp
  have hcond : ∀ u ∈ [[z.re, w.re]] ×ℂ [[z.im, w.im]],
      (u ≠ 1 → riemannZeta u ≠ 0) ∧ u ≠ 0 ∧ u ≠ -1 := by
    intro u hu
    rw [Complex.mem_reProdIm, hz_re, hw_re, hz_im, hw_im,
      Set.uIcc_of_le (by linarith : σ₁ ≤ σ₀), Set.uIcc_of_le (by linarith : -T ≤ T),
      Set.mem_Icc, Set.mem_Icc] at hu
    obtain ⟨⟨hre1, hre2⟩, him1, him2⟩ := hu
    refine ⟨?_, ?_, ?_⟩
    · intro hu1 hu0
      rcases lt_or_ge u.re 1 with hlt | hge
      · refine hzf u ?_ hlt (by rw [abs_le]; constructor <;> linarith) hu0
        rw [hσ₁def, hc₁def, hLdef] at hre1
        exact hre1
      · exact riemannZeta_ne_zero_of_one_le_re hge hu0
    · intro h0
      have := congrArg Complex.re h0
      simp at this
      linarith
    · intro hm1
      have := congrArg Complex.re hm1
      simp at this
      linarith
  have hshift := perron_F_vline_shift (by linarith : (0:ℝ) < x) z w
    (by rw [hz_re]; exact hσ₁lt) (by rw [hw_re]; exact hσ₀1)
    (by rw [hz_im]; linarith) (by rw [hw_im]; linarith) hcond
  rw [hz_re, hz_im, hw_re, hw_im] at hshift
  -- rearrange: RIGHT − 2π(x²/2) = LEFT + I·(bottom − top)
  have hrearr :
      (∫ t in (-T)..T,
          -logDeriv riemannZeta ((σ₀ : ℝ) + (t : ℝ) * Complex.I) *
              (x : ℂ) ^ (((σ₀ : ℝ) + (t : ℝ) * Complex.I) + 1) /
            (((σ₀ : ℝ) + (t : ℝ) * Complex.I) * (((σ₀ : ℝ) + (t : ℝ) * Complex.I) + 1)))
        - 2 * (Real.pi : ℂ) * ((x : ℂ) ^ (2 : ℂ) / 2)
      = (∫ t in (-T)..T,
          -logDeriv riemannZeta ((σ₁ : ℝ) + (t : ℝ) * Complex.I) *
              (x : ℂ) ^ (((σ₁ : ℝ) + (t : ℝ) * Complex.I) + 1) /
            (((σ₁ : ℝ) + (t : ℝ) * Complex.I) * (((σ₁ : ℝ) + (t : ℝ) * Complex.I) + 1)))
        + Complex.I * ((∫ t in σ₁..σ₀,
            -logDeriv riemannZeta ((t : ℝ) + ((-T : ℝ) : ℂ) * Complex.I) *
                (x : ℂ) ^ (((t : ℝ) + ((-T : ℝ) : ℂ) * Complex.I) + 1) /
              (((t : ℝ) + ((-T : ℝ) : ℂ) * Complex.I) *
                (((t : ℝ) + ((-T : ℝ) : ℂ) * Complex.I) + 1))) -
          (∫ t in σ₁..σ₀,
            -logDeriv riemannZeta ((t : ℝ) + ((T : ℝ) : ℂ) * Complex.I) *
                (x : ℂ) ^ (((t : ℝ) + ((T : ℝ) : ℂ) * Complex.I) + 1) /
              (((t : ℝ) + ((T : ℝ) : ℂ) * Complex.I) *
                (((t : ℝ) + ((T : ℝ) : ℂ) * Complex.I) + 1)))) := by
    linear_combination hshift
  rw [hrearr]
  -- bound the three pieces (typed in the L/σ₁-folded atoms; assignments are defeq-checked)
  have hleft' : ‖∫ t in (-T)..T,
      -logDeriv riemannZeta ((σ₁ : ℝ) + (t : ℝ) * Complex.I) *
          (x : ℂ) ^ (((σ₁ : ℝ) + (t : ℝ) * Complex.I) + 1) /
        (((σ₁ : ℝ) + (t : ℝ) * Complex.I) * (((σ₁ : ℝ) + (t : ℝ) * Complex.I) + 1))‖
      ≤ (192 * A + 4 * T₀ * C + 4 * T₀ / c₁) * (1 + L) * x ^ (σ₁ + 1) :=
    hlineb T hTT₁e x hx
  have htop := hconT T hTT₀ x hx σ₁ σ₀
    (by rw [hσ₁def, ← hLdef]) (by linarith) hσ₀2
  have hbot := hconB T hTT₀ x hx σ₁ σ₀
    (by rw [hσ₁def, ← hLdef]) (by linarith) hσ₀2
  have hσdiff : σ₀ - σ₁ ≤ 2 := by linarith
  have hkin : (0 : ℝ) ≤ A * L ^ 2 * x ^ (3 : ℝ) / T ^ 2 := by positivity
  have hconTle : ‖∫ t in σ₁..σ₀,
      -logDeriv riemannZeta ((t : ℝ) + ((T : ℝ) : ℂ) * Complex.I) *
          (x : ℂ) ^ (((t : ℝ) + ((T : ℝ) : ℂ) * Complex.I) + 1) /
        (((t : ℝ) + ((T : ℝ) : ℂ) * Complex.I) *
          (((t : ℝ) + ((T : ℝ) : ℂ) * Complex.I) + 1))‖
      ≤ A * L ^ 2 * x ^ (3 : ℝ) / T ^ 2 * 2 := by
    have h1 : ‖∫ t in σ₁..σ₀,
        -logDeriv riemannZeta ((t : ℝ) + ((T : ℝ) : ℂ) * Complex.I) *
            (x : ℂ) ^ (((t : ℝ) + ((T : ℝ) : ℂ) * Complex.I) + 1) /
          (((t : ℝ) + ((T : ℝ) : ℂ) * Complex.I) *
            (((t : ℝ) + ((T : ℝ) : ℂ) * Complex.I) + 1))‖
        ≤ A * L ^ 2 * x ^ (3 : ℝ) / T ^ 2 * (σ₀ - σ₁) := htop
    calc _ ≤ A * L ^ 2 * x ^ (3 : ℝ) / T ^ 2 * (σ₀ - σ₁) := h1
      _ ≤ A * L ^ 2 * x ^ (3 : ℝ) / T ^ 2 * 2 := by
          apply mul_le_mul_of_nonneg_left hσdiff hkin
  have hconBle : ‖∫ t in σ₁..σ₀,
      -logDeriv riemannZeta ((t : ℝ) + ((-T : ℝ) : ℂ) * Complex.I) *
          (x : ℂ) ^ (((t : ℝ) + ((-T : ℝ) : ℂ) * Complex.I) + 1) /
        (((t : ℝ) + ((-T : ℝ) : ℂ) * Complex.I) *
          (((t : ℝ) + ((-T : ℝ) : ℂ) * Complex.I) + 1))‖
      ≤ A * L ^ 2 * x ^ (3 : ℝ) / T ^ 2 * 2 := by
    have h1 : ‖∫ t in σ₁..σ₀,
        -logDeriv riemannZeta ((t : ℝ) + ((-T : ℝ) : ℂ) * Complex.I) *
            (x : ℂ) ^ (((t : ℝ) + ((-T : ℝ) : ℂ) * Complex.I) + 1) /
          (((t : ℝ) + ((-T : ℝ) : ℂ) * Complex.I) *
            (((t : ℝ) + ((-T : ℝ) : ℂ) * Complex.I) + 1))‖
        ≤ A * L ^ 2 * x ^ (3 : ℝ) / T ^ 2 * (σ₀ - σ₁) := hbot
    calc _ ≤ A * L ^ 2 * x ^ (3 : ℝ) / T ^ 2 * (σ₀ - σ₁) := h1
      _ ≤ A * L ^ 2 * x ^ (3 : ℝ) / T ^ 2 * 2 := by
          apply mul_le_mul_of_nonneg_left hσdiff hkin
  -- assemble
  apply le_trans (norm_add_le _ _)
  have hI2 : ‖Complex.I * ((∫ t in σ₁..σ₀,
      -logDeriv riemannZeta ((t : ℝ) + ((-T : ℝ) : ℂ) * Complex.I) *
          (x : ℂ) ^ (((t : ℝ) + ((-T : ℝ) : ℂ) * Complex.I) + 1) /
        (((t : ℝ) + ((-T : ℝ) : ℂ) * Complex.I) *
          (((t : ℝ) + ((-T : ℝ) : ℂ) * Complex.I) + 1))) -
      (∫ t in σ₁..σ₀,
      -logDeriv riemannZeta ((t : ℝ) + ((T : ℝ) : ℂ) * Complex.I) *
          (x : ℂ) ^ (((t : ℝ) + ((T : ℝ) : ℂ) * Complex.I) + 1) /
        (((t : ℝ) + ((T : ℝ) : ℂ) * Complex.I) *
          (((t : ℝ) + ((T : ℝ) : ℂ) * Complex.I) + 1))))‖
      ≤ A * L ^ 2 * x ^ (3 : ℝ) / T ^ 2 * 2 + A * L ^ 2 * x ^ (3 : ℝ) / T ^ 2 * 2 := by
    rw [norm_mul, Complex.norm_I, one_mul]
    exact le_trans (norm_sub_le _ _) (add_le_add hconBle hconTle)
  apply le_trans (add_le_add hleft' hI2)
  have h1 : (0 : ℝ) ≤ (1 + L) * x ^ (σ₁ + 1) := by positivity
  have h2 : (0 : ℝ) ≤ L ^ 2 * x ^ (3 : ℝ) / T ^ 2 := by positivity
  have hK1 : (0 : ℝ) ≤ 192 * A + 4 * T₀ * C + 4 * T₀ / c₁ := by
    have ha := mul_nonneg (by linarith : (0:ℝ) ≤ 4 * T₀) hC
    have hb := div_nonneg (by linarith : (0:ℝ) ≤ 4 * T₀) hc₁.le
    nlinarith [hA]
  have hexp : A * L ^ 2 * x ^ (3 : ℝ) / T ^ 2 * 2 + A * L ^ 2 * x ^ (3 : ℝ) / T ^ 2 * 2
      = 4 * A * (L ^ 2 * x ^ (3 : ℝ) / T ^ 2) := by ring
  rw [hexp]
  nlinarith [mul_nonneg hK1 h2, mul_nonneg (by linarith [hA] : (0:ℝ) ≤ 4 * A) h1]

#print axioms effective_truncated_perron
