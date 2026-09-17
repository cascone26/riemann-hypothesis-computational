import NormLSeriesLe
import Mathlib.NumberTheory.LSeries.Dirichlet

set_option maxHeartbeats 1000000

/-!
# §1591 — vertical-line bound: `‖L(Λ)(s)‖ ≤ ‖L(Λ)(Re s)‖` for `Re s > 1`

`‖term Λ s n‖` depends only on `Re s` (`norm_term_eq`), the terms at real `σ` are nonnegative reals, and the sum
of norms majorises the norm of the sum (§1421).  So on every vertical line `Re s = σ₀ > 1`, `-ζ'/ζ` is bounded by
its value at the real point — the uniform bound the `VerticalIntegrable` hypothesis of Mellin inversion needs.
EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction LSeries
open scoped LSeries.notation ComplexOrder

theorem norm_LSeries_vonMangoldt_le_line {s : ℂ} (hs : 1 < s.re) :
    ‖LSeries (↗Λ) s‖ ≤ ‖LSeries (↗Λ) ((s.re : ℝ) : ℂ)‖ := by
  have hσ : 1 < (((s.re : ℝ)) : ℂ).re := by rw [Complex.ofReal_re]; exact hs
  have hsum : Summable (fun n => ‖LSeries.term (↗Λ) s n‖) :=
    (ArithmeticFunction.LSeriesSummable_vonMangoldt hs).norm
  have hsumσ : LSeriesSummable (↗Λ) (((s.re : ℝ)) : ℂ) :=
    ArithmeticFunction.LSeriesSummable_vonMangoldt hσ
  have hterm : ∀ n, ‖LSeries.term (↗Λ) s n‖ = ‖LSeries.term (↗Λ) (((s.re : ℝ)) : ℂ) n‖ := by
    intro n
    rw [LSeries.norm_term_eq, LSeries.norm_term_eq, Complex.ofReal_re]
  have hnn : ∀ n, ‖LSeries.term (↗Λ) (((s.re : ℝ)) : ℂ) n‖
      = (LSeries.term (↗Λ) (((s.re : ℝ)) : ℂ) n).re := by
    intro n
    have hpos : (0 : ℂ) ≤ LSeries.term (↗Λ) (((s.re : ℝ)) : ℂ) n :=
      LSeries.term_nonneg (Complex.zero_le_real.mpr ArithmeticFunction.vonMangoldt_nonneg) s.re
    conv_rhs => rw [Complex.eq_coe_norm_of_nonneg hpos]
    rw [Complex.ofReal_re]
  calc ‖LSeries (↗Λ) s‖
      ≤ ∑' n, ‖LSeries.term (↗Λ) s n‖ := norm_LSeries_le hsum
    _ = ∑' n, ‖LSeries.term (↗Λ) (((s.re : ℝ)) : ℂ) n‖ := tsum_congr hterm
    _ = (LSeries (↗Λ) (((s.re : ℝ)) : ℂ)).re := by
        rw [tsum_congr hnn, (Complex.hasSum_re hsumσ.hasSum).tsum_eq]; rfl
    _ ≤ ‖LSeries (↗Λ) (((s.re : ℝ)) : ℂ)‖ := Complex.re_le_norm _

#print axioms norm_LSeries_vonMangoldt_le_line
