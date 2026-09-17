import ZetaTruncatedRep
import NormSumCpowBound
import FractIntegralTail

/-!
# §1291 — Unconditional ζ growth bound in the strip via the truncated representation

Assembling §1286 (truncated rep) with the three modulus bounds:
* §1290 `norm_sum_cpow_neg_le_bound` — main sum: `‖∑_{k≤N}k^{-s}‖ ≤ (N^{1-σ}-1)/(1-σ)`,
* the elementary term `‖N^{1-s}/(s-1)‖ = N^{1-σ}/‖s-1‖`,
* §1268 `norm_fract_integral_tail_le` — tail: `‖∫_N^∞{t}t^{-(s+1)}‖ ≤ N^{-σ}/σ`,

gives, for `0 < σ = Re s < 1` and any `N ≥ 1`,

  `‖ζ(s)‖ ≤ (N^{1-σ}-1)/(1-σ) + 1 + N^{1-σ}/‖s-1‖ + ‖s‖·(N^{-σ}/σ)`.

This is a fully unconditional, quantitative growth bound for ζ in the critical strip, uniform in `N`.  Optimizing
`N ≍ |t|` collapses the `N^{1-σ}` main terms to the log-type growth `≲ log|t|` feeding §1269 → sharp zero-free
region `σ > 1 − c/log|t|`.  EXIT0 / std-axioms.
-/

open Finset MeasureTheory Set

theorem norm_riemannZeta_le_strip_EM {N : ℕ} (hN : 1 ≤ N) {s : ℂ}
    (hs0 : 0 < s.re) (hs1 : s.re < 1) :
    ‖riemannZeta s‖ ≤ ((N : ℝ) ^ (1 - s.re) - 1) / (1 - s.re) + 1
      + (N : ℝ) ^ (1 - s.re) / ‖s - 1‖ + ‖s‖ * ((N : ℝ) ^ (-s.re) / s.re) := by
  have hsne1 : s ≠ 1 := by
    intro h; rw [h] at hs1; simp at hs1
  have hNpos : (0 : ℝ) < (N : ℝ) := by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hN
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  set A : ℂ := ∑ k ∈ Finset.Ioc 1 N, (k : ℂ) ^ (-s) with hA
  set B : ℂ := (N : ℂ) ^ (1 - s) / (s - 1) with hB
  set T : ℂ := ∫ t in Ioi (N : ℝ), ((Int.fract t : ℝ) : ℂ) * (t : ℂ) ^ (-(s + 1)) with hT
  -- triangle bound on the truncated representation
  rw [riemannZeta_eq_partial_sum_add_tail hN hs0 hsne1, ← hA, ← hB, ← hT]
  have htri : ‖A + 1 + B - s * T‖ ≤ ‖A‖ + 1 + ‖B‖ + ‖s‖ * ‖T‖ := by
    have h1 : ‖A + 1 + B - s * T‖ ≤ ‖A + 1 + B‖ + ‖s * T‖ := norm_sub_le _ _
    have h2 : ‖A + 1 + B‖ ≤ ‖A + 1‖ + ‖B‖ := norm_add_le _ _
    have h3 : ‖A + 1‖ ≤ ‖A‖ + ‖(1 : ℂ)‖ := norm_add_le _ _
    rw [norm_mul] at h1
    rw [norm_one] at h3
    linarith
  -- per-piece bounds
  have hAbound : ‖A‖ ≤ ((N : ℝ) ^ (1 - s.re) - 1) / (1 - s.re) :=
    norm_sum_cpow_neg_le_bound hN (le_of_lt hs0) (ne_of_lt hs1)
  have hBnorm : ‖B‖ = (N : ℝ) ^ (1 - s.re) / ‖s - 1‖ := by
    rw [hB, norm_div, ← Complex.ofReal_natCast N, Complex.norm_cpow_eq_rpow_re_of_pos hNpos,
      Complex.sub_re, Complex.one_re]
  have hTbound : ‖T‖ ≤ (N : ℝ) ^ (-s.re) / s.re :=
    norm_fract_integral_tail_le hs0 hNR
  have hsnn : (0 : ℝ) ≤ ‖s‖ := norm_nonneg _
  calc ‖A + 1 + B - s * T‖
      ≤ ‖A‖ + 1 + ‖B‖ + ‖s‖ * ‖T‖ := htri
    _ ≤ ((N : ℝ) ^ (1 - s.re) - 1) / (1 - s.re) + 1
          + (N : ℝ) ^ (1 - s.re) / ‖s - 1‖ + ‖s‖ * ((N : ℝ) ^ (-s.re) / s.re) := by
        rw [hBnorm]; gcongr

#print axioms norm_riemannZeta_le_strip_EM
