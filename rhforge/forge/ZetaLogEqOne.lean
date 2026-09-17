import ZetaTruncatedRep
import SumCpowNormBound
import SumHarmonicLe
import FractIntegralTail
import Term4SharpBound

/-!
# §1308 — ζ log bound on the line `σ = 1`: `‖ζ(1+it)‖ ≤ log N + 3 + 1/‖s-1‖`

The `σ = 1` case, which §1303/§1305 (needing `σ ≠ 1`) exclude but a Borel–Carathéodory disk centred at
`1 + it₀` passes through.  Assembled from §1286 (truncated rep, `s ≠ 1`), the harmonic main-sum bound
(§1287 ∘ §1307: `‖∑ k^{-s}‖ ≤ ∑ k⁻¹ ≤ log N`), and §1297 for the tail.  At `σ = 1` every `N`-power is `N⁰ = 1`:
`‖N^{1-s}/(s-1)‖ = 1/‖s-1‖` and `‖s‖·N^{-1} ≤ 2`, so `‖ζ(1+it)‖ ≤ log N + 3 + 1/‖s-1‖`.  EXIT0 / std-axioms.
-/

open Finset MeasureTheory Set

theorem norm_riemannZeta_le_log_eq_one {N : ℕ} (hN : 1 ≤ N) {s : ℂ}
    (hs1 : s.re = 1) (hsne : s ≠ 1) (hsim : |s.im| ≤ (N : ℝ)) :
    ‖riemannZeta s‖ ≤ Real.log N + 3 + 1 / ‖s - 1‖ := by
  have hs0 : 0 < s.re := by rw [hs1]; norm_num
  have hNpos : (0 : ℝ) < (N : ℝ) := by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hN
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  set A : ℂ := ∑ k ∈ Finset.Ioc 1 N, (k : ℂ) ^ (-s) with hA
  set B : ℂ := (N : ℂ) ^ (1 - s) / (s - 1) with hB
  set T : ℂ := ∫ t in Ioi (N : ℝ), ((Int.fract t : ℝ) : ℂ) * (t : ℂ) ^ (-(s + 1)) with hT
  rw [riemannZeta_eq_partial_sum_add_tail hN hs0 hsne, ← hA, ← hB, ← hT]
  have htri : ‖A + 1 + B - s * T‖ ≤ ‖A‖ + 1 + ‖B‖ + ‖s‖ * ‖T‖ := by
    have h1 := norm_sub_le (A + 1 + B) (s * T)
    have h2 := norm_add_le (A + 1) B
    have h3 := norm_add_le A (1 : ℂ)
    rw [norm_mul] at h1; rw [norm_one] at h3; linarith
  have hAbound : ‖A‖ ≤ Real.log N := by
    rw [hA]; refine le_trans norm_sum_cpow_neg_le ?_; rw [hs1]; exact sum_rpow_neg_one_le hN
  have hBnorm : ‖B‖ = 1 / ‖s - 1‖ := by
    rw [hB, norm_div, ← Complex.ofReal_natCast N, Complex.norm_cpow_eq_rpow_re_of_pos hNpos,
      Complex.sub_re, Complex.one_re, hs1, sub_self, Real.rpow_zero]
  have hpow0 : (N : ℝ) ^ (1 - s.re) / s.re = 1 := by
    rw [hs1, sub_self, Real.rpow_zero, div_one]
  have hTbound : ‖T‖ ≤ (N : ℝ) ^ (-s.re) / s.re := norm_fract_integral_tail_le hs0 hNR
  have hsn : (0 : ℝ) ≤ ‖s‖ := norm_nonneg _
  have hsT : ‖s‖ * ‖T‖ ≤ 2 := by
    have step1 : ‖s‖ * ‖T‖ ≤ ‖s‖ * ((N : ℝ) ^ (-s.re) / s.re) :=
      mul_le_mul_of_nonneg_left hTbound hsn
    have step2 : ‖s‖ * ((N : ℝ) ^ (-s.re) / s.re) ≤ 1 + (N : ℝ) ^ (1 - s.re) / s.re := by
      rw [← mul_div_assoc]; exact norm_mul_rpow_neg_div_le hs0 hNR hsim
    rw [hpow0] at step2
    linarith
  linarith [htri, hAbound, hBnorm.le, hsT]

#print axioms norm_riemannZeta_le_log_eq_one
