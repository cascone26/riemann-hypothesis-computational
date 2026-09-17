import NormZetaGeMoebiusTsum
import MoebiusTsumNormLe
import NormZetaOfRealLe

/-!
# §1443 ★★ — explicit ζ LOWER bound `‖ζ c‖ ≥ 1/(3 + 1/(Re c - 1))` for `1 < Re c`

Assembles the ζ lower-bound chain: §1422 (`‖ζ c‖ ≥ (∑'‖term ↗μ c n‖)⁻¹`) + §1442 (`∑' ≤ ‖ζ(σ)‖`) + §1441
(`‖ζ(σ)‖ ≤ 3+1/(σ-1)`), with `∑' ≥ ‖term ↗μ c 1‖ = 1 > 0` (the `n=1` term, `μ 1 = 1`).  Fully explicit `g₀`-free
lower bound in `σ = Re c` — the last input the BC scalar `hBexplicit` needs to be verified.  EXIT0 / std-axioms.
-/

open ArithmeticFunction
open scoped LSeries.notation Moebius

theorem norm_riemannZeta_ge_explicit {c : ℂ} (hc : 1 < c.re) :
    (3 + 1 / (c.re - 1))⁻¹ ≤ ‖riemannZeta c‖ := by
  have hlb := norm_riemannZeta_ge_moebius_tsum_inv hc
  have hub := moebius_tsum_norm_le hc
  have hζσ := norm_riemannZeta_ofReal_le hc
  have hμ : Summable (fun n => ‖LSeries.term ↗μ c n‖) :=
    summable_norm_iff.mpr (LSeriesSummable_moebius_iff.mpr hc)
  have hsum_le : ∑' n, ‖LSeries.term ↗μ c n‖ ≤ 3 + 1 / (c.re - 1) := le_trans hub hζσ
  have hterm1 : LSeries.term ↗μ c 1 = 1 := by
    rw [LSeries.term_of_ne_zero one_ne_zero]; simp
  have hsum_pos : 0 < ∑' n, ‖LSeries.term ↗μ c n‖ :=
    hμ.tsum_pos (fun n => norm_nonneg _) 1 (by rw [hterm1, norm_one]; exact one_pos)
  have hinv : (3 + 1 / (c.re - 1))⁻¹ ≤ (∑' n, ‖LSeries.term ↗μ c n‖)⁻¹ := by
    rw [inv_eq_one_div, inv_eq_one_div]; exact one_div_le_one_div_of_le hsum_pos hsum_le
  exact le_trans hinv hlb

#print axioms norm_riemannZeta_ge_explicit
