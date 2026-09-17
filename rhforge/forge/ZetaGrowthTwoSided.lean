import ZetaStripGrowth

/-!
# §1303 — Two-sided ζ growth bound (`σ ≠ 1`, both sides of the critical line's abscissa)

§1291 (`norm_riemannZeta_le_strip`) was stated for `0 < σ < 1`, but its ingredients only need `σ ≠ 1`:
`riemannZeta_eq_partial_sum_add_tail` (§1286) needs `s ≠ 1`; `norm_sum_cpow_neg_le_bound` (§1290) needs
`0 ≤ σ`, `σ ≠ 1`; `norm_fract_integral_tail_le` (§1268) needs `0 < σ`.  Relaxing `σ < 1` to `σ ≠ 1` gives the
same bound on BOTH sides of `σ = 1`:

  `‖ζ(s)‖ ≤ (N^{1-σ}-1)/(1-σ) + 1 + N^{1-σ}/‖s-1‖ + ‖s‖·(N^{-σ}/σ)`,   `0 < σ`, `σ ≠ 1`.

This is what a Borel–Carathéodory disk centred near `1 + it₀` (which straddles `σ = 1`) requires: a single ζ
bound valid across the whole disk.  EXIT0 / std-axioms.
-/

open Finset MeasureTheory Set

theorem norm_riemannZeta_le_offline {N : ℕ} (hN : 1 ≤ N) {s : ℂ}
    (hs0 : 0 < s.re) (hs1 : s.re ≠ 1) :
    ‖riemannZeta s‖ ≤ ((N : ℝ) ^ (1 - s.re) - 1) / (1 - s.re) + 1
      + (N : ℝ) ^ (1 - s.re) / ‖s - 1‖ + ‖s‖ * ((N : ℝ) ^ (-s.re) / s.re) := by
  have hsne1 : s ≠ 1 := fun h => hs1 (by rw [h]; exact Complex.one_re)
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
  have hAbound : ‖A‖ ≤ ((N : ℝ) ^ (1 - s.re) - 1) / (1 - s.re) :=
    norm_sum_cpow_neg_le_bound hN (le_of_lt hs0) hs1
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

#print axioms norm_riemannZeta_le_offline
