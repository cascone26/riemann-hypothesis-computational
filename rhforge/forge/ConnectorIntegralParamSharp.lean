import ConnectorPointwiseBound

set_option maxHeartbeats 2000000

open Complex Metric intervalIntegral

/-!
# §1702a — SHARP connector integral: x^{σ₀+1} in place of the crude x³: `‖∫ F‖ ≤ A·log²(T+2)·x³/T²·(σ₀-σ₁)` → 0

The final connector estimate: the Perron integrand `F = -ζ'/ζ(s)·x^{s+1}/(s(s+1))` on the
horizontal edge `Im = T` has `‖F‖ ≤ [A log²(T+2)]·x³/T²` (§1688 for the log-derivative factor,
`‖s‖,‖s+1‖ ≥ T` for the kinematic factor, `x^{σ+1} ≤ x³` for `x ≥ 1`, `σ ≤ 2`), so the edge
integral is `O(log²T/T²) → 0`. Combined with the §1668 shift-on-the-box identity this quantifies
the truncated Perron formula: main term `2π·x²/2` + explicitly-bounded connectors — the effective
smoothed-PNT inequality is now pure bookkeeping.
-/

theorem connector_integral_bound_param_sharp {c₁ A T₀ : ℝ}
    (hc₁ : 0 < c₁) (hA : 0 < A) (hT₀ : 6 ≤ T₀)
    (hpt : ∀ T : ℝ, T₀ ≤ T → ∀ σ : ℝ, 1 - c₁ / Real.log (T + 2) ≤ σ → σ ≤ 2 →
      ‖logDeriv riemannZeta ((σ : ℝ) + (T : ℝ) * Complex.I)‖
        ≤ A * (Real.log (T + 2)) ^ 2) :
    ∀ T : ℝ, T₀ ≤ T → ∀ x : ℝ, 1 ≤ x →
      ∀ σ₁ σ₀ : ℝ, 1 - c₁ / Real.log (T + 2) ≤ σ₁ → σ₁ ≤ σ₀ → σ₀ ≤ 2 →
      ‖∫ t in σ₁..σ₀,
          -logDeriv riemannZeta ((t : ℝ) + (T : ℝ) * Complex.I) *
              (x : ℂ) ^ (((t : ℝ) + (T : ℝ) * Complex.I) + 1) /
            (((t : ℝ) + (T : ℝ) * Complex.I) * (((t : ℝ) + (T : ℝ) * Complex.I) + 1))‖
        ≤ A * (Real.log (T + 2)) ^ 2 * x ^ (σ₀ + 1) / T ^ 2 * (σ₀ - σ₁) := by
  intro T hT x hx σ₁ σ₀ hσ₁lo hσle hσ₀hi
  have hT6 : (6 : ℝ) ≤ T := le_trans hT₀ hT
  have hTpos : (0 : ℝ) < T := by linarith
  have hL2 : (0 : ℝ) < Real.log (T + 2) := Real.log_pos (by linarith)
  have hbound : ∀ t ∈ Set.uIcc σ₁ σ₀,
      ‖-logDeriv riemannZeta ((t : ℝ) + (T : ℝ) * Complex.I) *
          (x : ℂ) ^ (((t : ℝ) + (T : ℝ) * Complex.I) + 1) /
        (((t : ℝ) + (T : ℝ) * Complex.I) * (((t : ℝ) + (T : ℝ) * Complex.I) + 1))‖
      ≤ A * (Real.log (T + 2)) ^ 2 * x ^ (σ₀ + 1) / T ^ 2 := by
    intro t ht
    rw [Set.uIcc_of_le hσle, Set.mem_Icc] at ht
    set s : ℂ := (t : ℝ) + (T : ℝ) * Complex.I with hsdef
    have hsre : s.re = t := by rw [hsdef]; simp
    have hsim : s.im = T := by rw [hsdef]; simp
    -- the log-derivative factor
    have hld : ‖logDeriv riemannZeta s‖ ≤ A * (Real.log (T + 2)) ^ 2 := by
      rw [hsdef]
      exact hpt T hT t (by linarith [ht.1]) (by linarith [ht.2])
    -- the power factor
    have hpow : ‖(x : ℂ) ^ (s + 1)‖ ≤ x ^ (σ₀ + 1) := by
      rw [Complex.norm_cpow_eq_rpow_re_of_pos (by linarith)]
      apply Real.rpow_le_rpow_of_exponent_le hx
      rw [Complex.add_re, hsre, Complex.one_re]
      linarith [ht.2]
    -- the denominator
    have hsnorm : T ≤ ‖s‖ := by
      calc T = |s.im| := by rw [hsim, abs_of_pos hTpos]
        _ ≤ ‖s‖ := abs_im_le_norm s
    have hs1norm : T ≤ ‖s + 1‖ := by
      calc T = |(s + 1).im| := by
            rw [show (s + 1).im = s.im by simp, hsim, abs_of_pos hTpos]
        _ ≤ ‖s + 1‖ := abs_im_le_norm _
    have hden : T ^ 2 ≤ ‖s * (s + 1)‖ := by
      rw [norm_mul]
      calc T ^ 2 = T * T := sq T
        _ ≤ ‖s‖ * ‖s + 1‖ := by
            apply mul_le_mul hsnorm hs1norm hTpos.le (norm_nonneg _)
    -- assemble
    rw [norm_div, norm_mul, norm_neg]
    have hnum : ‖logDeriv riemannZeta s‖ * ‖(x : ℂ) ^ (s + 1)‖
        ≤ A * (Real.log (T + 2)) ^ 2 * x ^ (σ₀ + 1) := by
      apply mul_le_mul hld hpow (norm_nonneg _) (by positivity)
    have hdenpos : (0 : ℝ) < ‖s * (s + 1)‖ := lt_of_lt_of_le (by positivity) hden
    have hstep1 : ‖logDeriv riemannZeta s‖ * ‖(x : ℂ) ^ (s + 1)‖ / ‖s * (s + 1)‖
        ≤ (A * (Real.log (T + 2)) ^ 2 * x ^ (σ₀ + 1)) / ‖s * (s + 1)‖ := by
      gcongr
    have hstep2 : (A * (Real.log (T + 2)) ^ 2 * x ^ (σ₀ + 1)) / ‖s * (s + 1)‖
        ≤ A * (Real.log (T + 2)) ^ 2 * x ^ (σ₀ + 1) / T ^ 2 := by
      gcongr
    exact hstep1.trans hstep2
  calc ‖∫ t in σ₁..σ₀, -logDeriv riemannZeta ((t : ℝ) + (T : ℝ) * Complex.I) *
          (x : ℂ) ^ (((t : ℝ) + (T : ℝ) * Complex.I) + 1) /
        (((t : ℝ) + (T : ℝ) * Complex.I) * (((t : ℝ) + (T : ℝ) * Complex.I) + 1))‖
      ≤ A * (Real.log (T + 2)) ^ 2 * x ^ (σ₀ + 1) / T ^ 2 * |σ₀ - σ₁| :=
        intervalIntegral.norm_integral_le_of_norm_le_const
          (fun t ht => hbound t (Set.uIoc_subset_uIcc ht))
    _ = A * (Real.log (T + 2)) ^ 2 * x ^ (σ₀ + 1) / T ^ 2 * (σ₀ - σ₁) := by
        rw [abs_of_nonneg (by linarith)]

#print axioms connector_integral_bound_param_sharp
