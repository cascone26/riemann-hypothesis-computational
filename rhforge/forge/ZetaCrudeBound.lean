import ZetaGrowthTwoSided
import SumCpowNormBound
import SumRpowNegLeCard

/-!
# §1439 — crude ζ bound UNIFORM across `Re s = 1`

The offline bound §1303 degenerates at `σ = 1` (the `(N^{1-σ}-1)/(1-σ)` term).  Replacing that sharp partial-sum
bound with the crude uniform `‖∑_{k∈Ioc 1 N} k^{-s}‖ ≤ ∑ k^{-σ} ≤ N` (§`norm_sum_cpow_neg_le` + §1438) removes the
`σ ≠ 1` requirement — the resulting bound holds for ALL `0 < Re s`, `s ≠ 1`:

  `‖ζ(s)‖ ≤ N + 1 + N^{1-σ}/‖s-1‖ + ‖s‖·(N^{-σ}/σ)`.

This is exactly what a Borel–Carathéodory sphere straddling `Re s = 1` needs: one explicit bound valid at every
point.  EXIT0 / std-axioms.
-/

open Finset MeasureTheory Set

theorem norm_riemannZeta_le_crude {N : ℕ} (hN : 1 ≤ N) {s : ℂ} (hs0 : 0 < s.re) (hsne1 : s ≠ 1) :
    ‖riemannZeta s‖ ≤ (N : ℝ) + 1
      + (N : ℝ) ^ (1 - s.re) / ‖s - 1‖ + ‖s‖ * ((N : ℝ) ^ (-s.re) / s.re) := by
  have hNpos : (0 : ℝ) < (N : ℝ) := by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hN
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  set A : ℂ := ∑ k ∈ Finset.Ioc 1 N, (k : ℂ) ^ (-s) with hA
  set B : ℂ := (N : ℂ) ^ (1 - s) / (s - 1) with hB
  set T : ℂ := ∫ t in Ioi (N : ℝ), ((Int.fract t : ℝ) : ℂ) * (t : ℂ) ^ (-(s + 1)) with hT
  rw [riemannZeta_eq_partial_sum_add_tail hN hs0 hsne1, ← hA, ← hB, ← hT]
  have htri : ‖A + 1 + B - s * T‖ ≤ ‖A‖ + 1 + ‖B‖ + ‖s‖ * ‖T‖ := by
    have h1 : ‖A + 1 + B - s * T‖ ≤ ‖A + 1 + B‖ + ‖s * T‖ := norm_sub_le _ _
    have h2 : ‖A + 1 + B‖ ≤ ‖A + 1‖ + ‖B‖ := norm_add_le _ _
    have h3 : ‖A + 1‖ ≤ ‖A‖ + ‖(1 : ℂ)‖ := norm_add_le _ _
    rw [norm_mul] at h1
    rw [norm_one] at h3
    linarith
  have hAbound : ‖A‖ ≤ (N : ℝ) :=
    le_trans norm_sum_cpow_neg_le (sum_rpow_neg_le_card (le_of_lt hs0))
  have hBnorm : ‖B‖ = (N : ℝ) ^ (1 - s.re) / ‖s - 1‖ := by
    rw [hB, norm_div, ← Complex.ofReal_natCast N, Complex.norm_cpow_eq_rpow_re_of_pos hNpos,
      Complex.sub_re, Complex.one_re]
  have hTbound : ‖T‖ ≤ (N : ℝ) ^ (-s.re) / s.re :=
    norm_fract_integral_tail_le hs0 hNR
  have hsnn : (0 : ℝ) ≤ ‖s‖ := norm_nonneg _
  calc ‖A + 1 + B - s * T‖
      ≤ ‖A‖ + 1 + ‖B‖ + ‖s‖ * ‖T‖ := htri
    _ ≤ (N : ℝ) + 1 + (N : ℝ) ^ (1 - s.re) / ‖s - 1‖ + ‖s‖ * ((N : ℝ) ^ (-s.re) / s.re) := by
        rw [hBnorm]; gcongr

#print axioms norm_riemannZeta_le_crude
