import HpoleRiemannZetaBox
import LogDerivZetaVonMangoldt
import NormLSeriesVonMangoldtLine

set_option maxHeartbeats 1000000

/-!
# §1703b — THE POLE BOUND ON `‖L(Λ)(σ₀)‖`: `∃ A ≥ 0, ‖L(Λ)(σ₀)‖ ≤ 1/(σ₀-1) + A` on `1 < σ₀ ≤ 2`

At a REAL point `σ₀` the von Mangoldt series has nonnegative real terms, so its norm IS its real
part; and its real part is `Re(-ζ'/ζ(σ₀))`, which the box pole bound (§1572) caps by
`1/(σ₀-1) + A`.  Feeds the §1703 optimization: at `σ₀ = 1 + 1/log x` this is `log x + A`.
-/

open Complex ArithmeticFunction LSeries
open scoped LSeries.notation ComplexOrder

theorem norm_LSeries_vonMangoldt_pole_bound :
    ∃ A : ℝ, 0 ≤ A ∧ ∀ σ₀ : ℝ, 1 < σ₀ → σ₀ ≤ 2 →
      ‖LSeries (↗Λ) ((σ₀ : ℝ) : ℂ)‖ ≤ 1 / (σ₀ - 1) + A := by
  obtain ⟨A, hA0, hA⟩ := hpole_riemannZeta_box
  refine ⟨A, hA0, fun σ₀ h1 h2 => ?_⟩
  have hσre : (((σ₀ : ℝ)) : ℂ).re = σ₀ := Complex.ofReal_re _
  have hσ : 1 < (((σ₀ : ℝ)) : ℂ).re := by rw [hσre]; exact h1
  have hsumσ : LSeriesSummable (↗Λ) (((σ₀ : ℝ)) : ℂ) :=
    ArithmeticFunction.LSeriesSummable_vonMangoldt hσ
  have hnn : ∀ n, ‖LSeries.term (↗Λ) (((σ₀ : ℝ)) : ℂ) n‖
      = (LSeries.term (↗Λ) (((σ₀ : ℝ)) : ℂ) n).re := by
    intro n
    have hpos : (0 : ℂ) ≤ LSeries.term (↗Λ) (((σ₀ : ℝ)) : ℂ) n :=
      LSeries.term_nonneg (Complex.zero_le_real.mpr ArithmeticFunction.vonMangoldt_nonneg) σ₀
    conv_rhs => rw [Complex.eq_coe_norm_of_nonneg hpos]
    rw [Complex.ofReal_re]
  have hLre : ‖LSeries (↗Λ) ((σ₀ : ℝ) : ℂ)‖ ≤ (LSeries (↗Λ) ((σ₀ : ℝ) : ℂ)).re := by
    calc ‖LSeries (↗Λ) ((σ₀ : ℝ) : ℂ)‖
        ≤ ∑' n, ‖LSeries.term (↗Λ) ((σ₀ : ℝ) : ℂ) n‖ := norm_LSeries_le hsumσ.norm
      _ = (LSeries (↗Λ) ((σ₀ : ℝ) : ℂ)).re := by
          rw [tsum_congr hnn, (Complex.hasSum_re hsumσ.hasSum).tsum_eq]; rfl
  have hLeq : LSeries (↗Λ) ((σ₀ : ℝ) : ℂ) = -logDeriv riemannZeta ((σ₀ : ℝ) : ℂ) := by
    have h := logDeriv_riemannZeta_eq_neg_vonMangoldt hσ
    linear_combination h
  have hbox := hA ((σ₀ : ℝ) : ℂ) hσ (by rw [hσre]; exact h2)
    (by rw [Complex.ofReal_im]; simp)
  rw [hσre] at hbox
  calc ‖LSeries (↗Λ) ((σ₀ : ℝ) : ℂ)‖
      ≤ (LSeries (↗Λ) ((σ₀ : ℝ) : ℂ)).re := hLre
    _ = (-logDeriv riemannZeta ((σ₀ : ℝ) : ℂ)).re := by rw [hLeq]
    _ ≤ 1 / (σ₀ - 1) + A := hbox

#print axioms norm_LSeries_vonMangoldt_pole_bound
