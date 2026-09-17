import ConnectorPointwiseBound
import LogDerivZetaConj
import LogDerivZetaSmallT

set_option maxHeartbeats 2000000

open Complex Metric

/-!
# §1692 ★★ — the LEFT-LINE POINTWISE BOUND: `ζ'/ζ` on the shifted line `Re = 1 - c₁/log(T+2)`

The full-height estimate on the zero-free vertical line, three regimes glued:
- `t ≥ T₁`: the §1688 connector bound APPLIED AT HEIGHT `t` (valid since the shifted line sits
  right of every height-`t` zero-free line: `σ₁(T) ≥ 1 - c₁/log(t+2)` for `t ≤ T`), giving
  `A·log²(t+2)` — the `t`-LOCAL log factor that makes the error integral converge;
- `t ≤ -T₁`: free by the §1690 Schwarz reflection;
- `|t| ≤ T₁`: the §1691 pole-split compact bound `C + log(T+2)/c₁`.
This is the integrand bound for the left-line error integral (§1693): the `x^{σ₁+1}`-weighted
integral is `O(x^{σ₁+1}·(1 + log(T+2)))` — the classical error-term shape.
-/

theorem left_line_pointwise_bound :
    ∃ c₁ A T₁ Tsp C : ℝ, 0 < c₁ ∧ 0 < A ∧ 6 ≤ T₁ ∧ 6 ≤ Tsp ∧ Tsp ≤ T₁ ∧ 0 ≤ C ∧
      ∀ T : ℝ, T₁ ≤ T → ∀ t : ℝ, |t| ≤ T →
      (Tsp ≤ |t| →
        ‖logDeriv riemannZeta (((1 - c₁ / Real.log (T + 2) : ℝ) : ℂ) + (t : ℝ) * Complex.I)‖
          ≤ A * (Real.log (|t| + 2)) ^ 2) ∧
      (|t| ≤ Tsp →
        ‖logDeriv riemannZeta (((1 - c₁ / Real.log (T + 2) : ℝ) : ℂ) + (t : ℝ) * Complex.I)‖
          ≤ C + Real.log (T + 2) / c₁) := by
  obtain ⟨c₁, A, T₀, hc₁, hA, hT₀, hconn⟩ := connector_pointwise_bound
  obtain ⟨σ₂, C, hσ₂half, hσ₂lt, hC, hsmall⟩ := logDeriv_zeta_small_t_bound T₀ (by linarith)
  set T₁ : ℝ := max T₀ (Real.exp (c₁ / (1 - σ₂)) + 6) with hT₁def
  have hT₁T₀ : T₀ ≤ T₁ := le_max_left _ _
  have hT₁exp : Real.exp (c₁ / (1 - σ₂)) + 6 ≤ T₁ := le_max_right _ _
  have hT₁6 : (6 : ℝ) ≤ T₁ := le_trans hT₀ hT₁T₀
  refine ⟨c₁, A, T₁, T₀, C, hc₁, hA, hT₁6, hT₀, hT₁T₀, hC, ?_⟩
  intro T hT t ht
  have hT6 : (6 : ℝ) ≤ T := le_trans hT₁6 hT
  have hLpos : 0 < Real.log (T + 2) := Real.log_pos (by linarith)
  set σ₁ : ℝ := 1 - c₁ / Real.log (T + 2) with hσ₁def
  have hσ₁lt : σ₁ < 1 := by
    rw [hσ₁def]
    have : 0 < c₁ / Real.log (T + 2) := by positivity
    linarith
  have hσ₁gtσ₂ : σ₂ < σ₁ := by
    have hlogbig : c₁ / (1 - σ₂) < Real.log (T + 2) := by
      rw [← Real.log_exp (c₁ / (1 - σ₂))]
      apply Real.log_lt_log (Real.exp_pos _)
      have := Real.exp_pos (c₁ / (1 - σ₂))
      linarith
    have h1σ₂ : 0 < 1 - σ₂ := by linarith
    have hkey : c₁ / Real.log (T + 2) < 1 - σ₂ := by
      rw [div_lt_iff₀ hLpos]
      calc c₁ = c₁ / (1 - σ₂) * (1 - σ₂) := by field_simp
        _ < Real.log (T + 2) * (1 - σ₂) := mul_lt_mul_of_pos_right hlogbig h1σ₂
        _ = (1 - σ₂) * Real.log (T + 2) := by ring
    rw [hσ₁def]
    linarith
  constructor
  · -- large |t| via the height-|t| connector bound (+ reflection for t < 0)
    intro htbig
    have htabs6 : (6 : ℝ) ≤ |t| := le_trans hT₀ htbig
    have hσ₁range : 1 - c₁ / Real.log (|t| + 2) ≤ σ₁ := by
      rw [hσ₁def]
      have hlt : Real.log (|t| + 2) ≤ Real.log (T + 2) :=
        Real.log_le_log (by positivity) (by linarith)
      have hLtpos : 0 < Real.log (|t| + 2) := Real.log_pos (by linarith)
      have : c₁ / Real.log (T + 2) ≤ c₁ / Real.log (|t| + 2) :=
        div_le_div_of_nonneg_left hc₁.le hLtpos hlt
      linarith
    rcases le_or_gt 0 t with htpos | htneg
    · have habs : |t| = t := abs_of_nonneg htpos
      rw [habs] at htabs6 hσ₁range ⊢
      exact hconn t (by rw [← habs]; exact htbig) σ₁ hσ₁range (hσ₁lt.le.trans (by norm_num))
    · have habs : |t| = -t := abs_of_neg htneg
      have hconj : (((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I)
          = (starRingEnd ℂ) (((σ₁ : ℝ) : ℂ) + ((-t : ℝ) : ℂ) * Complex.I) := by
        rw [map_add, map_mul]
        simp [Complex.conj_I]
      have hne1 : (((σ₁ : ℝ) : ℂ) + ((-t : ℝ) : ℂ) * Complex.I) ≠ 1 := by
        intro h
        have := congrArg Complex.im h
        simp at this
        linarith
      rw [hconj, norm_logDeriv_riemannZeta_conj hne1]
      rw [habs] at htabs6 hσ₁range ⊢
      exact hconn (-t) (by rw [← habs]; exact htbig) σ₁ hσ₁range
        (hσ₁lt.le.trans (by norm_num))
  · -- small |t| via the pole-split compact bound
    intro htsmall
    have hs1 : (((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I) ≠ 1 := by
      intro h
      have := congrArg Complex.re h
      simp at this
      linarith
    have hsre : (((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I).re = σ₁ := by simp
    have hsim : (((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I).im = t := by simp
    have hb := hsmall (((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I)
      (by rw [hsre]; exact hσ₁gtσ₂) (by rw [hsre]; linarith)
      (by rw [hsim]; exact htsmall) hs1
    -- the pole term on the shifted line
    have hpole : 1 / ‖(((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I) - 1‖ ≤ Real.log (T + 2) / c₁ := by
      have hre : ((((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I) - 1).re = σ₁ - 1 := by simp
      have hnorm : 1 - σ₁ ≤ ‖(((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I) - 1‖ := by
        calc 1 - σ₁ = |σ₁ - 1| := by rw [abs_sub_comm, abs_of_nonneg (by linarith)]
          _ = |((((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I) - 1).re| := by rw [hre]
          _ ≤ ‖(((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I) - 1‖ := Complex.abs_re_le_norm _
      have h1σ₁ : 0 < 1 - σ₁ := by linarith
      have hnormpos : 0 < ‖(((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I) - 1‖ :=
        lt_of_lt_of_le h1σ₁ hnorm
      have hstep : 1 / ‖(((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I) - 1‖ ≤ 1 / (1 - σ₁) :=
        one_div_le_one_div_of_le h1σ₁ hnorm
      have heq : 1 / (1 - σ₁) = Real.log (T + 2) / c₁ := by
        rw [hσ₁def]
        rw [show (1 : ℝ) - (1 - c₁ / Real.log (T + 2)) = c₁ / Real.log (T + 2) by ring,
          one_div_div]
      linarith
    linarith

#print axioms left_line_pointwise_bound
