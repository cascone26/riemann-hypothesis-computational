import ZetaCrudeBound
import SumHarmonicLe
import ZetaImagAxisExact

/-!
# S2977 — ζ on Re s = 1 and the convexity endpoint |ζ(it)| ≪ |t|^{1/2} log|t| (lane 1)

* `norm_riemannZeta_le_re_one` — Re s = 1, s ≠ 1, N ≥ 1: ‖ζ s‖ ≤ log N + 1 + 1/‖s−1‖ + ‖s‖/N
  (partial sum Σ_{1<k≤N} k^{-1} ≤ log N + Euler–Maclaurin tail, valid ON the line Re s = 1).
* `norm_zeta_one_sub_mul_I_le` — t ≥ 2: ‖ζ(1 − it)‖ ≤ log t + 5   (N = ⌊t⌋).
* `norm_zeta_mul_I_sq_le_log` — t ≥ 2: ‖ζ(it)‖² ≤ (t/2π)(log t + 5)²   (with S2975).
The left edge of the convexity strip is now explicit; Phragmén–Lindelöf across 0 ≤ Re s ≤ 1 is the next brick.
-/

open Complex Finset MeasureTheory Set

namespace ZetaImagAxisBound

theorem norm_riemannZeta_le_re_one {N : ℕ} (hN : 1 ≤ N) {s : ℂ} (hs1 : s.re = 1) (hsne : s ≠ 1) :
    ‖riemannZeta s‖ ≤ Real.log N + 1 + 1 / ‖s - 1‖ + ‖s‖ / N := by
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  have hNpos : (0 : ℝ) < (N : ℝ) := by linarith
  have hs0 : 0 < s.re := by rw [hs1]; norm_num
  set A : ℂ := ∑ k ∈ Finset.Ioc 1 N, (k : ℂ) ^ (-s) with hA
  set B : ℂ := (N : ℂ) ^ (1 - s) / (s - 1) with hB
  set T : ℂ := ∫ t in Ioi (N : ℝ), ((Int.fract t : ℝ) : ℂ) * (t : ℂ) ^ (-(s + 1)) with hT
  rw [riemannZeta_eq_partial_sum_add_tail hN hs0 hsne, ← hA, ← hB, ← hT]
  have htri : ‖A + 1 + B - s * T‖ ≤ ‖A‖ + 1 + ‖B‖ + ‖s‖ * ‖T‖ := by
    have h1 : ‖A + 1 + B - s * T‖ ≤ ‖A + 1 + B‖ + ‖s * T‖ := norm_sub_le _ _
    have h2 : ‖A + 1 + B‖ ≤ ‖A + 1‖ + ‖B‖ := norm_add_le _ _
    have h3 : ‖A + 1‖ ≤ ‖A‖ + ‖(1 : ℂ)‖ := norm_add_le _ _
    rw [norm_mul] at h1
    rw [norm_one] at h3
    linarith
  have hAb : ‖A‖ ≤ Real.log N := by
    refine le_trans norm_sum_cpow_neg_le ?_
    rw [hs1]
    exact sum_rpow_neg_one_le hN
  have hBn : ‖B‖ = 1 / ‖s - 1‖ := by
    rw [hB, norm_div, ← Complex.ofReal_natCast N, Complex.norm_cpow_eq_rpow_re_of_pos hNpos,
      Complex.sub_re, Complex.one_re, hs1]
    simp
  have hTb : ‖T‖ ≤ (N : ℝ) ^ (-s.re) / s.re := norm_fract_integral_tail_le hs0 hNR
  rw [hs1, Real.rpow_neg_one, div_one] at hTb
  have hsT : ‖s‖ * ‖T‖ ≤ ‖s‖ / N := by
    rw [div_eq_mul_inv]; exact mul_le_mul_of_nonneg_left hTb (norm_nonneg _)
  linarith

theorem norm_zeta_one_sub_mul_I_le {t : ℝ} (ht : 2 ≤ t) :
    ‖riemannZeta (1 - (t : ℂ) * I)‖ ≤ Real.log t + 5 := by
  have hN : 1 ≤ ⌊t⌋₊ := Nat.le_floor (by norm_num; linarith)
  have hs1 : (1 - (t : ℂ) * I).re = 1 := by simp
  have hsne : (1 - (t : ℂ) * I) ≠ 1 := by
    intro h; have := congrArg Complex.im h; simp at this; linarith
  have h := norm_riemannZeta_le_re_one hN hs1 hsne
  have hNx : (⌊t⌋₊ : ℝ) ≤ t := Nat.floor_le (by linarith)
  have hNlt : t < (⌊t⌋₊ : ℝ) + 1 := Nat.lt_floor_add_one t
  have hNpos : (0 : ℝ) < (⌊t⌋₊ : ℝ) := by linarith
  have hlog : Real.log (⌊t⌋₊ : ℝ) ≤ Real.log t := Real.log_le_log hNpos hNx
  have hs1n : ‖(1 - (t : ℂ) * I) - 1‖ = t := by
    rw [show (1 - (t : ℂ) * I) - 1 = -((t : ℂ) * I) by ring, norm_neg, norm_mul, Complex.norm_I,
      mul_one, Complex.norm_real, Real.norm_of_nonneg (by linarith)]
  have hsn : ‖(1 - (t : ℂ) * I)‖ ≤ 1 + t := by
    calc ‖(1 - (t : ℂ) * I)‖ ≤ ‖(1 : ℂ)‖ + ‖(t : ℂ) * I‖ := norm_sub_le _ _
      _ = 1 + t := by
        rw [norm_one, norm_mul, Complex.norm_I, mul_one, Complex.norm_real,
          Real.norm_of_nonneg (by linarith)]
  rw [hs1n] at h
  have h1t : 1 / t ≤ 1 := by rw [div_le_one (by linarith)]; linarith
  have h2 : ‖(1 - (t : ℂ) * I)‖ / (⌊t⌋₊ : ℝ) ≤ 3 := by
    rw [div_le_iff₀ hNpos]; linarith
  linarith

theorem norm_zeta_mul_I_sq_le_log {t : ℝ} (ht : 2 ≤ t) :
    ‖riemannZeta ((t : ℂ) * I)‖ ^ 2 ≤ t / (2 * Real.pi) * (Real.log t + 5) ^ 2 := by
  have h := ZetaImagAxis.norm_zeta_mul_I_sq_le t (by linarith)
  have hb := norm_zeta_one_sub_mul_I_le ht
  have hlog : 0 ≤ Real.log t := Real.log_nonneg (by linarith)
  have hsq : ‖riemannZeta (1 - (t : ℂ) * I)‖ ^ 2 ≤ (Real.log t + 5) ^ 2 :=
    pow_le_pow_left₀ (norm_nonneg _) hb 2
  have hc : 0 ≤ t / (2 * Real.pi) := by positivity
  exact le_trans h (mul_le_mul_of_nonneg_left hsq hc)

end ZetaImagAxisBound

#print axioms ZetaImagAxisBound.norm_riemannZeta_le_re_one
#print axioms ZetaImagAxisBound.norm_zeta_one_sub_mul_I_le
#print axioms ZetaImagAxisBound.norm_zeta_mul_I_sq_le_log
