import NormLSeriesVonMangoldtLine
import Mathlib.NumberTheory.LSeries.Dirichlet

set_option maxHeartbeats 1000000

/-!
# §1628 — monotonicity of `‖L(Λ)(σ)‖` in `σ`: `σ₁ ≤ σ₂ ⟹ ‖L(Λ)(σ₂)‖ ≤ ‖L(Λ)(σ₁)‖`  (`σ₁ > 1`)

The von Mangoldt L-series at a real point `σ > 1` is a sum of nonnegative reals `∑ Λ(n)/n^σ`, so its norm equals its
real part `∑ (term)ᵣₑ`, and each term is decreasing in `σ` (`n^{-σ}` decreasing, `Λ ≥ 0`).  Hence `‖L(Λ)(σ)‖` is
decreasing in `σ`.  This gives the uniform bound `‖L(Λ)(σ')‖ ≤ ‖L(Λ)(σ₁)‖` over a strip `σ' ∈ [σ₁, σ₀]`, needed to
show the horizontal edges of the shift rectangle integrate to `0`.  EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction LSeries
open scoped LSeries.notation ComplexOrder

theorem norm_LSeries_vonMangoldt_mono {σ₁ σ₂ : ℝ} (hσ₁ : 1 < σ₁) (hle : σ₁ ≤ σ₂) :
    ‖LSeries (↗Λ) ((σ₂ : ℝ) : ℂ)‖ ≤ ‖LSeries (↗Λ) ((σ₁ : ℝ) : ℂ)‖ := by
  have hσ₂ : 1 < σ₂ := lt_of_lt_of_le hσ₁ hle
  have hs₁ : LSeriesSummable (↗Λ) ((σ₁ : ℝ) : ℂ) :=
    ArithmeticFunction.LSeriesSummable_vonMangoldt (by rwa [Complex.ofReal_re])
  have hs₂ : LSeriesSummable (↗Λ) ((σ₂ : ℝ) : ℂ) :=
    ArithmeticFunction.LSeriesSummable_vonMangoldt (by rwa [Complex.ofReal_re])
  -- termwise nonnegativity of the real part = the norm
  have hterm_re : ∀ (σ : ℝ) (n : ℕ),
      (LSeries.term (↗Λ) ((σ : ℝ) : ℂ) n).re = ‖LSeries.term (↗Λ) ((σ : ℝ) : ℂ) n‖ := by
    intro σ n
    have hpos : (0 : ℂ) ≤ LSeries.term (↗Λ) ((σ : ℝ) : ℂ) n :=
      LSeries.term_nonneg (Complex.zero_le_real.mpr ArithmeticFunction.vonMangoldt_nonneg) σ
    conv_lhs => rw [Complex.eq_coe_norm_of_nonneg hpos]
    rw [Complex.ofReal_re]
  -- `‖L(Λ)(σ)‖ = ∑' (term σ n).re`
  have hnorm : ∀ (σ : ℝ), LSeriesSummable (↗Λ) ((σ : ℝ) : ℂ) →
      ‖LSeries (↗Λ) ((σ : ℝ) : ℂ)‖ = ∑' n, (LSeries.term (↗Λ) ((σ : ℝ) : ℂ) n).re := by
    intro σ hsum
    have hpos : (0 : ℂ) ≤ LSeries (↗Λ) ((σ : ℝ) : ℂ) := by
      rw [LSeries]
      exact tsum_nonneg (fun n =>
        LSeries.term_nonneg (Complex.zero_le_real.mpr ArithmeticFunction.vonMangoldt_nonneg) σ)
    have hre : ‖LSeries (↗Λ) ((σ : ℝ) : ℂ)‖ = (LSeries (↗Λ) ((σ : ℝ) : ℂ)).re := by
      conv_rhs => rw [Complex.eq_coe_norm_of_nonneg hpos]
      rw [Complex.ofReal_re]
    rw [hre]
    exact ((Complex.hasSum_re hsum.hasSum).tsum_eq).symm
  rw [hnorm σ₂ hs₂, hnorm σ₁ hs₁]
  refine (Complex.hasSum_re hs₂.hasSum).summable.tsum_le_tsum (fun n => ?_)
    (Complex.hasSum_re hs₁.hasSum).summable
  rw [hterm_re σ₂ n, hterm_re σ₁ n, LSeries.norm_term_eq, LSeries.norm_term_eq,
    Complex.ofReal_re, Complex.ofReal_re]
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp
  · rw [if_neg hn.ne', if_neg hn.ne']
    apply div_le_div_of_nonneg_left
    · exact norm_nonneg _
    · exact Real.rpow_pos_of_pos (Nat.cast_pos.mpr hn) _
    · exact Real.rpow_le_rpow_of_exponent_le (Nat.one_le_cast.mpr hn) hle

#print axioms norm_LSeries_vonMangoldt_mono
