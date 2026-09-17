import MoebiusTsumNormLe

set_option maxHeartbeats 1000000

/-!
# §1537 — `∑' ‖term (↗χ·↗μ) c n‖ ≤ ‖ζ(Re c)‖` for `1 < Re c`

L-analogue of ζ §1442 — but lighter: rather than re-running the `1`-series real-part machinery, bound the twist
TERMWISE by the plain Möbius series (`‖(↗χ·↗μ) n‖ = ‖χ(n)‖·‖μ(n)‖ ≤ ‖μ(n)‖` via `χ.norm_le_one`), then delegate to
§1442 (`∑' ‖term ↗μ c n‖ ≤ ‖ζ σ‖`).  With §1536 and §1441 this yields `‖L(c,χ)‖ ≥ 1/(3+1/(σ-1))` — the same bound
as ζ, because `|χ(n)| ≤ 1` bounds the twist by the plain Möbius majorant.  EXIT0 / std-axioms.
-/

open ArithmeticFunction Complex DirichletCharacter
open scoped LSeries.notation ArithmeticFunction.Moebius

theorem twisted_moebius_tsum_norm_le {N : ℕ} (χ : DirichletCharacter ℂ N) {c : ℂ} (hc : 1 < c.re) :
    ∑' n, ‖LSeries.term (↗χ * ↗μ) c n‖ ≤ ‖riemannZeta ((c.re : ℝ) : ℂ)‖ := by
  have hχμ := summable_norm_iff.mpr (LSeriesSummable_mul χ (LSeriesSummable_moebius_iff.mpr hc))
  have hμ := summable_norm_iff.mpr
    (LSeriesSummable_moebius_iff (s := c).mpr hc)
  have hcmp : ∀ n, ‖LSeries.term (↗χ * ↗μ) c n‖ ≤ ‖LSeries.term (↗μ) c n‖ := by
    intro n
    refine LSeries.norm_term_le c ?_
    rw [Pi.mul_apply, norm_mul]
    exact mul_le_of_le_one_left (norm_nonneg _) (χ.norm_le_one _)
  calc ∑' n, ‖LSeries.term (↗χ * ↗μ) c n‖
      ≤ ∑' n, ‖LSeries.term (↗μ) c n‖ := hasSum_le hcmp hχμ.hasSum hμ.hasSum
    _ ≤ ‖riemannZeta ((c.re : ℝ) : ℂ)‖ := moebius_tsum_norm_le hc

#print axioms twisted_moebius_tsum_norm_le
