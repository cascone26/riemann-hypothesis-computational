import Mathlib.NumberTheory.LSeries.Dirichlet
import Mathlib.Analysis.PSeries

set_option maxHeartbeats 1000000

open Complex LSeries ArithmeticFunction

open scoped LSeries.notation ArithmeticFunction.Moebius

/-!
# §1678 — THE ζ LOWER BOUND: `‖ζ(s)‖ ≥ 1/(Σ n^{-σ})` for `Re s = σ > 1`

The Möbius-series lower bound: `ζ(s)·L(μ)(s) = 1` and `‖L(μ)(s)‖ ≤ Σ|μ(n)|/n^σ ≤ Σ 1/n^σ`, so
`‖ζ(s)‖ ≥ (Σ 1/n^σ)⁻¹`. At the Landau center `z₀ = 1 + 1/log T + iT` this gives
`‖ζ(z₀)‖ ≥ 1/O(log T)` — the `‖f c‖`-input of §1677's master bound (`log(M/‖ζ(z₀)‖) = O(log T)`).
-/

theorem inv_tsum_le_norm_riemannZeta {s : ℂ} (hs : 1 < s.re) :
    (∑' n : ℕ, 1 / (n : ℝ) ^ s.re)⁻¹ ≤ ‖riemannZeta s‖ := by
  -- the inverse identity
  have hmul : riemannZeta s * L ↗μ s = 1 := by
    rw [← LSeries_zeta_eq_riemannZeta hs]
    exact LSeries_zeta_mul_Lseries_moebius hs
  -- summability facts
  have hμsum : LSeriesSummable ↗μ s := ArithmeticFunction.LSeriesSummable_moebius_iff.mpr hs
  have hpsum : Summable (fun n : ℕ => 1 / (n : ℝ) ^ s.re) :=
    Real.summable_one_div_nat_rpow.mpr hs
  have hppos : 0 < ∑' n : ℕ, 1 / (n : ℝ) ^ s.re := by
    have h1le : 1 / ((1 : ℕ) : ℝ) ^ s.re ≤ ∑' n : ℕ, 1 / (n : ℝ) ^ s.re :=
      hpsum.le_tsum 1 (fun j _ => by positivity)
    simp only [Nat.cast_one, Real.one_rpow, div_one] at h1le
    linarith
  -- ‖L μ s‖ ≤ Σ 1/n^σ
  have hLμ : ‖L ↗μ s‖ ≤ ∑' n : ℕ, 1 / (n : ℝ) ^ s.re := by
    calc ‖L ↗μ s‖ ≤ ∑' n : ℕ, ‖term ↗μ s n‖ := norm_tsum_le_tsum_norm hμsum.norm
      _ ≤ ∑' n : ℕ, 1 / (n : ℝ) ^ s.re := by
          apply hμsum.norm.tsum_le_tsum _ hpsum
          intro n
          rw [LSeries.norm_term_eq]
          rcases eq_or_ne n 0 with rfl | hn
          · simp
            positivity
          · rw [if_neg hn]
            rw [div_le_div_iff_of_pos_right]
            · have : ‖(↗μ n : ℂ)‖ = |(μ n : ℝ)| := by
                rw [show ((↗μ) n : ℂ) = ((μ n : ℝ) : ℂ) by push_cast; rfl,
                  Complex.norm_real, Real.norm_eq_abs]
              rw [this]
              have := ArithmeticFunction.abs_moebius_le_one (n := n)
              calc |(μ n : ℝ)| = ((|μ n| : ℤ) : ℝ) := by push_cast; ring
                _ ≤ 1 := by exact_mod_cast this
            · have : 0 < (n : ℝ) := by positivity
              positivity
  -- combine
  have hζnorm : ‖riemannZeta s‖ * ‖L ↗μ s‖ = 1 := by
    rw [← norm_mul, hmul, norm_one]
  have hLμpos : 0 < ‖L ↗μ s‖ := by
    by_contra hle
    push_neg at hle
    have : ‖L ↗μ s‖ = 0 := le_antisymm hle (norm_nonneg _)
    rw [this, mul_zero] at hζnorm
    norm_num at hζnorm
  rw [inv_le_iff_one_le_mul₀ hppos]
  calc 1 = ‖riemannZeta s‖ * ‖L ↗μ s‖ := hζnorm.symm
    _ ≤ ‖riemannZeta s‖ * (∑' n : ℕ, 1 / (n : ℝ) ^ s.re) := by
        gcongr

#print axioms inv_tsum_le_norm_riemannZeta
