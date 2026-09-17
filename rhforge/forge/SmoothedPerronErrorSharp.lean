import PsiPerronBridge
import PerronTailBound
import EffectiveTruncatedPerronSharp

set_option maxHeartbeats 4000000

open Complex MeasureTheory Chebyshev

open scoped LSeries.notation ArithmeticFunction

/-!
# §1702d — the free-T smoothed PNT, SHARP (x^{σ₀+1} everywhere)

For all `T ≥ T₂`, `x ≥ 1`, `σ₀ ∈ (1, 2]`:

  `‖2π·Ψ(x) − 2π·x²/2‖ ≤ K((1+log(T+2))·x^{σ₁(T)+1} + log²(T+2)·x³/T² + ‖L(Λ)(σ₀)‖·x^{σ₀+1}/T)`

— the SMOOTHED PRIME NUMBER THEOREM WITH ERROR TERM at every finite truncation height, fully
machine-checked from tonight's chain: `2πΨ = ∫_ℝ F` (§1699 Mellin bridge), the tail `∫_ℝ − ∫_T`
(§1700), and the truncated Perron `∫_T = 2π·x²/2 + error` (§1698, itself the winding arc + the
Landau ladder + the zero-free region). Taking `σ₀ = 1 + 1/log x`, `T = e^{√log x}` — a pure
arithmetic optimization pass — gives `Ψ(x) = x²/2 + O(x²e^{-c√log x})`, i.e. de la Vallée
Poussin's smoothed PNT. Every constant here is explicit or a banked existential.
-/

theorem smoothed_perron_error_free_T_sharp :
    ∃ c₁ K T₂ : ℝ, 0 < c₁ ∧ 0 < K ∧ 6 ≤ T₂ ∧ ∀ T : ℝ, T₂ ≤ T → ∀ x : ℝ, 1 ≤ x →
      ∀ σ₀ : ℝ, 1 < σ₀ → σ₀ ≤ 2 →
      ‖(2 * Real.pi : ℝ) • psi1C x - 2 * (Real.pi : ℂ) * ((x : ℂ) ^ (2 : ℂ) / 2)‖
        ≤ K * ((1 + Real.log (T + 2)) * x ^ (1 - c₁ / Real.log (T + 2) + 1)
            + (Real.log (T + 2)) ^ 2 * x ^ (σ₀ + 1) / T ^ 2
            + ‖LSeries (↗Λ) ((σ₀ : ℝ) : ℂ)‖ * x ^ (σ₀ + 1) / T) := by
  obtain ⟨c₁, K, T₂, hc₁, hK, hT₂, hperron⟩ := effective_truncated_perron_sharp
  refine ⟨c₁, K + 2, T₂, hc₁, by linarith, hT₂, ?_⟩
  intro T hT x hx σ₀ hσ₀1 hσ₀2
  have hT6 : (6 : ℝ) ≤ T := le_trans hT₂ hT
  have hxpos : (0 : ℝ) < x := by linarith
  -- the three ingredients
  have hbridge := two_pi_psi1C_eq_perron_integral hσ₀1 hxpos
  have htail := perron_tail_bound hx hσ₀1 (T := T) (by linarith)
  have htrunc := hperron T hT x hx σ₀ hσ₀1 hσ₀2
  -- triangle: 2πΨ − main = (∫_ℝ − ∫_T) + (∫_T − main)
  have hdecomp : (2 * Real.pi : ℝ) • psi1C x - 2 * (Real.pi : ℂ) * ((x : ℂ) ^ (2 : ℂ) / 2)
      = ((∫ t : ℝ,
          -logDeriv riemannZeta ((σ₀ : ℂ) + (t : ℝ) * I) *
              (x : ℂ) ^ (((σ₀ : ℂ) + (t : ℝ) * I) + 1) /
            (((σ₀ : ℂ) + (t : ℝ) * I) * (((σ₀ : ℂ) + (t : ℝ) * I) + 1)))
        - ∫ t in (-T)..T,
          -logDeriv riemannZeta ((σ₀ : ℂ) + (t : ℝ) * I) *
              (x : ℂ) ^ (((σ₀ : ℂ) + (t : ℝ) * I) + 1) /
            (((σ₀ : ℂ) + (t : ℝ) * I) * (((σ₀ : ℂ) + (t : ℝ) * I) + 1)))
      + ((∫ t in (-T)..T,
          -logDeriv riemannZeta ((σ₀ : ℂ) + (t : ℝ) * I) *
              (x : ℂ) ^ (((σ₀ : ℂ) + (t : ℝ) * I) + 1) /
            (((σ₀ : ℂ) + (t : ℝ) * I) * (((σ₀ : ℂ) + (t : ℝ) * I) + 1)))
        - 2 * (Real.pi : ℂ) * ((x : ℂ) ^ (2 : ℂ) / 2)) := by
    rw [hbridge]
    ring
  rw [hdecomp]
  have hLpos : (0 : ℝ) < Real.log (T + 2) := Real.log_pos (by linarith)
  calc ‖_ + _‖ ≤ _ + _ := norm_add_le _ _
    _ ≤ (2 * ‖LSeries (↗Λ) ((σ₀ : ℝ) : ℂ)‖ * x ^ (σ₀ + 1) / T)
        + K * ((1 + Real.log (T + 2)) * x ^ (1 - c₁ / Real.log (T + 2) + 1)
            + (Real.log (T + 2)) ^ 2 * x ^ (σ₀ + 1) / T ^ 2) :=
        add_le_add htail htrunc
    _ ≤ (K + 2) * ((1 + Real.log (T + 2)) * x ^ (1 - c₁ / Real.log (T + 2) + 1)
            + (Real.log (T + 2)) ^ 2 * x ^ (σ₀ + 1) / T ^ 2
            + ‖LSeries (↗Λ) ((σ₀ : ℝ) : ℂ)‖ * x ^ (σ₀ + 1) / T) := by
        have h1 : (0 : ℝ) ≤ (1 + Real.log (T + 2)) * x ^ (1 - c₁ / Real.log (T + 2) + 1) := by
          positivity
        have h2 : (0 : ℝ) ≤ (Real.log (T + 2)) ^ 2 * x ^ (σ₀ + 1) / T ^ 2 := by
          positivity
        have h3 : (0 : ℝ) ≤ ‖LSeries (↗Λ) ((σ₀ : ℝ) : ℂ)‖ * x ^ (σ₀ + 1) / T := by
          have := norm_nonneg (LSeries (↗Λ) ((σ₀ : ℝ) : ℂ))
          have hxp : (0 : ℝ) ≤ x ^ (σ₀ + 1) := (Real.rpow_pos_of_pos hxpos _).le
          positivity
        have he : 2 * ‖LSeries (↗Λ) ((σ₀ : ℝ) : ℂ)‖ * x ^ (σ₀ + 1) / T
            = 2 * (‖LSeries (↗Λ) ((σ₀ : ℝ) : ℂ)‖ * x ^ (σ₀ + 1) / T) := by ring
        rw [he]
        nlinarith [mul_nonneg hK.le h3, mul_nonneg hK.le h1, mul_nonneg hK.le h2, h1, h2, h3]

#print axioms smoothed_perron_error_free_T_sharp
