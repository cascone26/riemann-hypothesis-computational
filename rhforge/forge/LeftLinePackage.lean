import ConnectorPointwiseBound
import LogDerivZetaConj
import LogDerivZetaSmallT
import LandauZeroDistance

set_option maxHeartbeats 2000000

open Complex Metric

/-!
# §1695 ★★ — the LEFT-LINE PACKAGE: one coherent constant for bounds + zero-freeness

Supersedes §1692's over-specialized form. The subtlety: the pointwise bounds (§1688/§1691) and the
zero-distance property (§1686) each produce their own existential constant, and separately-obtained
existentials cannot be proven equal. Fix: obtain BOTH and take `c₁ := min` — every property survives
for the smaller constant, since shrinking `c₁` moves the band `Re ≥ 1 - c₁/log(T+2)` RIGHT, deeper
into each regime's validity region. One package, flexible `σ`: for `σ ∈ [1 - c₁/log(T+2), 1)`,
`|t| ≤ T`: the `t`-local `A·log²(|t|+2)` bound (`|t| ≥ Tsp`), the pole-split `C + 1/(1-σ)` bound
(`|t| ≤ Tsp`), and ζ-freeness of the whole band below `Re = 1`.
-/

theorem left_line_package :
    ∃ c₁ A T₁ Tsp C : ℝ, 0 < c₁ ∧ 0 < A ∧ 6 ≤ T₁ ∧ 6 ≤ Tsp ∧ Tsp ≤ T₁ ∧ 0 ≤ C ∧
      ∀ T : ℝ, T₁ ≤ T →
      (∀ σ : ℝ, 1 - c₁ / Real.log (T + 2) ≤ σ → σ < 1 → ∀ t : ℝ, |t| ≤ T →
        (Tsp ≤ |t| →
          ‖logDeriv riemannZeta (((σ : ℝ) : ℂ) + (t : ℝ) * Complex.I)‖
            ≤ A * (Real.log (|t| + 2)) ^ 2) ∧
        (|t| ≤ Tsp →
          ‖logDeriv riemannZeta (((σ : ℝ) : ℂ) + (t : ℝ) * Complex.I)‖
            ≤ C + 1 / (1 - σ))) ∧
      (∀ z : ℂ, 1 - c₁ / Real.log (T + 2) ≤ z.re → z.re < 1 → |z.im| ≤ T + 1 →
        riemannZeta z ≠ 0) := by
  obtain ⟨cA, A, T₀, hcA, hA, hT₀, hconn⟩ := connector_pointwise_bound
  obtain ⟨cB, hcB, hdist⟩ := landau_ball_zero_distance
  obtain ⟨σ₂, C, hσ₂half, hσ₂lt, hC, hsmall⟩ := logDeriv_zeta_small_t_bound T₀ (by linarith)
  set c₁ : ℝ := min cA cB with hc₁def
  have hc₁ : 0 < c₁ := lt_min hcA hcB
  have hc₁A : c₁ ≤ cA := min_le_left _ _
  have hc₁B : c₁ ≤ cB := min_le_right _ _
  set T₁ : ℝ := max T₀ (Real.exp (c₁ / (1 - σ₂)) + 6) with hT₁def
  have hT₁T₀ : T₀ ≤ T₁ := le_max_left _ _
  have hT₁exp : Real.exp (c₁ / (1 - σ₂)) + 6 ≤ T₁ := le_max_right _ _
  have hT₁6 : (6 : ℝ) ≤ T₁ := le_trans hT₀ hT₁T₀
  refine ⟨c₁, A, T₁, T₀, C, hc₁, hA, hT₁6, hT₀, hT₁T₀, hC, ?_⟩
  intro T hT
  have hT6 : (6 : ℝ) ≤ T := le_trans hT₁6 hT
  have hLpos : 0 < Real.log (T + 2) := Real.log_pos (by linarith)
  constructor
  · intro σ hσlo hσhi t ht
    have hσσ₂ : σ₂ < σ := by
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
      linarith
    constructor
    · -- large |t|: the height-|t| connector bound (reflection for t < 0)
      intro htbig
      have htabs6 : (6 : ℝ) ≤ |t| := le_trans hT₀ htbig
      have hσrange : 1 - cA / Real.log (|t| + 2) ≤ σ := by
        have hlt : Real.log (|t| + 2) ≤ Real.log (T + 2) :=
          Real.log_le_log (by positivity) (by linarith)
        have hLtpos : 0 < Real.log (|t| + 2) := Real.log_pos (by linarith)
        have h1 : c₁ / Real.log (T + 2) ≤ cA / Real.log (|t| + 2) := by
          calc c₁ / Real.log (T + 2) ≤ c₁ / Real.log (|t| + 2) :=
                div_le_div_of_nonneg_left hc₁.le hLtpos hlt
            _ ≤ cA / Real.log (|t| + 2) := by gcongr
        linarith
      rcases le_or_gt 0 t with htpos | htneg
      · have habs : |t| = t := abs_of_nonneg htpos
        rw [habs] at htabs6 hσrange ⊢
        exact hconn t (by rw [← habs]; exact htbig) σ hσrange (by linarith)
      · have habs : |t| = -t := abs_of_neg htneg
        have hconj : (((σ : ℝ) : ℂ) + (t : ℝ) * Complex.I)
            = (starRingEnd ℂ) (((σ : ℝ) : ℂ) + ((-t : ℝ) : ℂ) * Complex.I) := by
          rw [map_add, map_mul]
          simp [Complex.conj_I]
        have hne1 : (((σ : ℝ) : ℂ) + ((-t : ℝ) : ℂ) * Complex.I) ≠ 1 := by
          intro h
          have := congrArg Complex.im h
          simp at this
          linarith
        rw [hconj, norm_logDeriv_riemannZeta_conj hne1]
        rw [habs] at htabs6 hσrange ⊢
        exact hconn (-t) (by rw [← habs]; exact htbig) σ hσrange (by linarith)
    · -- small |t|: the pole-split bound
      intro htsmall
      have hs1 : (((σ : ℝ) : ℂ) + (t : ℝ) * Complex.I) ≠ 1 := by
        intro h
        have := congrArg Complex.re h
        simp at this
        linarith
      have hsre : (((σ : ℝ) : ℂ) + (t : ℝ) * Complex.I).re = σ := by simp
      have hsim : (((σ : ℝ) : ℂ) + (t : ℝ) * Complex.I).im = t := by simp
      have hb := hsmall (((σ : ℝ) : ℂ) + (t : ℝ) * Complex.I)
        (by rw [hsre]; exact hσσ₂) (by rw [hsre]; linarith)
        (by rw [hsim]; exact htsmall) hs1
      have hpole : 1 / ‖(((σ : ℝ) : ℂ) + (t : ℝ) * Complex.I) - 1‖ ≤ 1 / (1 - σ) := by
        have hre : ((((σ : ℝ) : ℂ) + (t : ℝ) * Complex.I) - 1).re = σ - 1 := by simp
        have hnorm : 1 - σ ≤ ‖(((σ : ℝ) : ℂ) + (t : ℝ) * Complex.I) - 1‖ := by
          calc 1 - σ = |σ - 1| := by rw [abs_sub_comm, abs_of_nonneg (by linarith)]
            _ = |((((σ : ℝ) : ℂ) + (t : ℝ) * Complex.I) - 1).re| := by rw [hre]
            _ ≤ ‖(((σ : ℝ) : ℂ) + (t : ℝ) * Complex.I) - 1‖ := Complex.abs_re_le_norm _
        exact one_div_le_one_div_of_le (by linarith) hnorm
      linarith
  · -- zero-freeness of the band
    intro z hzre hzre1 hzim
    intro hz0
    have hzreB : 1 - cB / Real.log (T + 2) ≤ z.re := by
      have : c₁ / Real.log (T + 2) ≤ cB / Real.log (T + 2) := by gcongr
      linarith
    have h := hdist T (by linarith) z z hzreB hzim hz0
    rw [sub_self, norm_zero] at h
    have : 0 < cB / Real.log (T + 2) := by positivity
    linarith

#print axioms left_line_package
