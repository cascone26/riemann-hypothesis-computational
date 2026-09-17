import Mathlib.NumberTheory.LSeries.Dirichlet

/-!
# §1442 — `∑' ‖term ↗μ c n‖ ≤ ‖ζ(Re c)‖` for `1 < Re c`

Termwise `‖term ↗μ c n‖ ≤ ‖term 1 c n‖` (`|μ n| ≤ 1`, same `s = c`), then `‖term 1 c n‖ = ‖term 1 (Re c) n‖`
(depends only on the real part).  The majorant sums to `(ζ σ).re ≤ ‖ζ σ‖` (the `1`-series terms are nonnegative
reals).  With §1422 and §1441 this gives the explicit ζ lower bound the BC scalar needs.  EXIT0 / std-axioms.
-/

open ArithmeticFunction Complex
open scoped LSeries.notation Moebius ComplexOrder

theorem moebius_tsum_norm_le {c : ℂ} (hc : 1 < c.re) :
    ∑' n, ‖LSeries.term ↗μ c n‖ ≤ ‖riemannZeta ((c.re : ℝ) : ℂ)‖ := by
  have hσc : 1 < (((c.re : ℝ)) : ℂ).re := by rw [Complex.ofReal_re]; exact hc
  have hμ : Summable (fun n => ‖LSeries.term ↗μ c n‖) :=
    summable_norm_iff.mpr (LSeriesSummable_moebius_iff.mpr hc)
  have h1sum : LSeriesSummable (1 : ℕ → ℂ) (((c.re : ℝ)) : ℂ) := LSeriesSummable_one_iff.mpr hσc
  have h1 : Summable (fun n => ‖LSeries.term (1 : ℕ → ℂ) (((c.re : ℝ)) : ℂ) n‖) :=
    summable_norm_iff.mpr h1sum
  have hcmp : ∀ n, ‖(↗μ) n‖ ≤ ‖(1 : ℕ → ℂ) n‖ := by
    intro n
    rw [Pi.one_apply, norm_one, Complex.norm_intCast]
    exact_mod_cast abs_moebius_le_one
  have hterm : ∀ n, ‖LSeries.term ↗μ c n‖ ≤ ‖LSeries.term (1 : ℕ → ℂ) (((c.re : ℝ)) : ℂ) n‖ := by
    intro n
    have hst : ‖LSeries.term (1 : ℕ → ℂ) c n‖ = ‖LSeries.term (1 : ℕ → ℂ) (((c.re : ℝ)) : ℂ) n‖ := by
      rw [LSeries.norm_term_eq, LSeries.norm_term_eq, Complex.ofReal_re]
    rw [← hst]; exact LSeries.norm_term_le c (hcmp n)
  have hnn : ∀ n, ‖LSeries.term (1 : ℕ → ℂ) (((c.re : ℝ)) : ℂ) n‖
      = (LSeries.term (1 : ℕ → ℂ) (((c.re : ℝ)) : ℂ) n).re := by
    intro n
    have hpos : (0 : ℂ) ≤ LSeries.term (1 : ℕ → ℂ) (((c.re : ℝ)) : ℂ) n :=
      LSeries.term_nonneg (by rw [Pi.one_apply]; exact zero_le_one) c.re
    conv_rhs => rw [Complex.eq_coe_norm_of_nonneg hpos]
    rw [Complex.ofReal_re]
  have heq : ∑' n, ‖LSeries.term (1 : ℕ → ℂ) (((c.re : ℝ)) : ℂ) n‖
      = (riemannZeta (((c.re : ℝ)) : ℂ)).re := by
    rw [tsum_congr hnn, (Complex.hasSum_re h1sum.hasSum).tsum_eq]
    congr 1
    exact LSeries_one_eq_riemannZeta hσc
  calc ∑' n, ‖LSeries.term ↗μ c n‖
      ≤ ∑' n, ‖LSeries.term (1 : ℕ → ℂ) (((c.re : ℝ)) : ℂ) n‖ :=
        hasSum_le hterm hμ.hasSum h1.hasSum
    _ = (riemannZeta (((c.re : ℝ)) : ℂ)).re := heq
    _ ≤ ‖riemannZeta (((c.re : ℝ)) : ℂ)‖ := Complex.re_le_norm _

#print axioms moebius_tsum_norm_le
