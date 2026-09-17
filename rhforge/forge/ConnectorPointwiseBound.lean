import ZetaLogDerivSqBound
import NormLSeriesVonMangoldtMono
import NormLSeriesVonMangoldtLine

set_option maxHeartbeats 2000000

open Complex Metric ArithmeticFunction LSeries

open scoped LSeries.notation

/-!
# §1688 ★★ — the CONNECTOR POINTWISE BOUND: `‖ζ'/ζ(σ+iT)‖ ≤ A·log²(T+2)` for `σ ∈ [σ₁(T), 2]`

The full horizontal-edge bound, both regimes: near the 1-line (`σ ≤ 1+1/Lg+1/64`) via the §1687
Landau-ball bound; away (`σ > 1+1/64`) via the Dirichlet series `ζ'/ζ = -L(Λ)` with §1591's
line-domination and §1628's monotonicity, giving the CONSTANT `‖L(Λ)(1+1/64)‖`. `T₀ = e^{64(1+c₁)}`
makes the two regimes overlap. This is the exact per-point input for the §1666/§1668 connector
integrals: with `F`'s kinematic factor `x^{s+1}/(s(s+1))` (which contributes `1/T²`, §1627-pattern),
the connectors are `O(x^{σ₀+1}·log²T/T²) → 0`.
-/

theorem connector_pointwise_bound :
    ∃ c₁ A T₀ : ℝ, 0 < c₁ ∧ 0 < A ∧ 6 ≤ T₀ ∧ ∀ T : ℝ, T₀ ≤ T → ∀ σ : ℝ,
      1 - c₁ / Real.log (T + 2) ≤ σ → σ ≤ 2 →
      ‖logDeriv riemannZeta ((σ : ℝ) + (T : ℝ) * Complex.I)‖
        ≤ A * (Real.log (T + 2)) ^ 2 := by
  obtain ⟨c₁, A₁, hc₁, hA₁, hball⟩ := norm_logDeriv_zeta_sq_bound
  set C₂ : ℝ := ‖LSeries (↗Λ) (((1 + 1 / 64 : ℝ) : ℝ) : ℂ)‖
    with hC₂def
  refine ⟨c₁, A₁ + C₂, Real.exp (64 * (1 + c₁)) + 6, hc₁, ?_, ?_, ?_⟩
  · have : 0 ≤ C₂ := norm_nonneg _
    linarith
  · have := Real.exp_pos (64 * (1 + c₁))
    linarith
  intro T hT σ hσlo hσhi
  set Lg : ℝ := Real.log (T + 2) with hLgdef
  have hT6 : (6 : ℝ) ≤ T := by
    have := Real.exp_pos (64 * (1 + c₁))
    linarith
  have hLpos : 0 < Lg := by
    rw [hLgdef]
    exact Real.log_pos (by linarith)
  have hL2 : 2 ≤ Lg := by
    rw [hLgdef, Real.le_log_iff_exp_le (by linarith)]
    have he := Real.exp_one_lt_d9
    calc Real.exp 2 = Real.exp 1 * Real.exp 1 := by rw [← Real.exp_add]; norm_num
      _ ≤ 2.7182818286 * 2.7182818286 :=
          mul_le_mul he.le he.le (Real.exp_pos 1).le (by norm_num)
      _ ≤ T + 2 := by nlinarith
  have hLbig : 64 * (1 + c₁) ≤ Lg := by
    rw [hLgdef, Real.le_log_iff_exp_le (by linarith)]
    linarith
  have hAL : 0 ≤ (Real.log (T + 2)) ^ 2 := sq_nonneg _
  by_cases hcase : σ ≤ 1 + 1 / Lg + 1 / 64
  · -- Landau-ball regime
    have hzdist : ‖((σ : ℝ) : ℂ) + (T : ℝ) * Complex.I
        - (((1 + 1 / Lg : ℝ) : ℂ) + (T : ℝ) * Complex.I)‖ ≤ 1 / 64 := by
      have heq : ((σ : ℝ) : ℂ) + (T : ℝ) * Complex.I
          - (((1 + 1 / Lg : ℝ) : ℂ) + (T : ℝ) * Complex.I)
          = ((σ - (1 + 1 / Lg) : ℝ) : ℂ) := by
        push_cast
        ring
      rw [heq, Complex.norm_real, Real.norm_eq_abs, abs_le]
      constructor
      · -- σ - (1+1/Lg) ≥ -1/64: from σ ≥ 1 - c₁/Lg and (1+c₁)/Lg ≤ 1/64
        have hdiv : (1 + c₁) / Lg ≤ 1 / 64 := by
          rw [div_le_iff₀ hLpos]
          linarith
        have : c₁ / Lg + 1 / Lg ≤ 1 / 64 := by
          rw [show c₁ / Lg + 1 / Lg = (1 + c₁) / Lg by ring]
          exact hdiv
        have h1Lg : 0 < 1 / Lg := by positivity
        linarith
      · linarith
    have hzre : 1 - c₁ / Lg ≤ (((σ : ℝ) : ℂ) + (T : ℝ) * Complex.I).re := by
      simpa using hσlo
    have := hball T hT6 (((σ : ℝ) : ℂ) + (T : ℝ) * Complex.I) hzdist hzre
    calc ‖logDeriv riemannZeta (((σ : ℝ) : ℂ) + (T : ℝ) * Complex.I)‖
        ≤ A₁ * (Real.log (T + 2)) ^ 2 := this
      _ ≤ (A₁ + C₂) * (Real.log (T + 2)) ^ 2 := by
          have : 0 ≤ C₂ := norm_nonneg _
          nlinarith
  · -- Dirichlet-series regime: σ > 1 + 1/64
    push_neg at hcase
    have hσgt : 1 + 1 / 64 ≤ σ := by
      have h1Lg : 0 < 1 / Lg := by positivity
      linarith
    have hsre : (((σ : ℝ) : ℂ) + (T : ℝ) * Complex.I).re = σ := by simp
    have hsre1 : 1 < (((σ : ℝ) : ℂ) + (T : ℝ) * Complex.I).re := by
      rw [hsre]
      linarith
    have hld : logDeriv riemannZeta (((σ : ℝ) : ℂ) + (T : ℝ) * Complex.I)
        = -LSeries (↗Λ) (((σ : ℝ) : ℂ) + (T : ℝ) * Complex.I) := by
      rw [logDeriv_riemannZeta_eq_neg_vonMangoldt hsre1]
    rw [hld, norm_neg]
    have h1 := norm_LSeries_vonMangoldt_le_line hsre1
    rw [hsre] at h1
    have h2 := norm_LSeries_vonMangoldt_mono (σ₁ := 1 + 1 / 64) (σ₂ := σ)
      (by norm_num) hσgt
    calc ‖LSeries (↗Λ) (((σ : ℝ) : ℂ) + (T : ℝ) * Complex.I)‖
        ≤ ‖LSeries (↗Λ) ((σ : ℝ) : ℂ)‖ := h1
      _ ≤ C₂ := h2
      _ ≤ C₂ * (Real.log (T + 2)) ^ 2 := by
          have hLsq : 1 ≤ (Real.log (T + 2)) ^ 2 := by
            rw [← hLgdef]
            nlinarith
          nlinarith [norm_nonneg (LSeries (↗Λ) (((1 + 1 / 64 : ℝ) : ℝ) : ℂ))]
      _ ≤ (A₁ + C₂) * (Real.log (T + 2)) ^ 2 := by
          nlinarith

#print axioms connector_pointwise_bound
