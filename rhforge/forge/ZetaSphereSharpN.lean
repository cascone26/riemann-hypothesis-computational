import ZetaCrudeBound

/-!
# S2959 — sharp N-term ζ bound uniform over a sphere (lane-1 structural lever)

The 5-8-4-1 pipeline (S2945–S2957) feeds Borel–Carathéodory with the N = 1 sphere bound
(`norm_riemannZeta_le_sphere_N1`), whose `‖c‖/(c.re-R)` term is ≍ |t|, so log M ≍ log|t|.
This brick replaces the crude partial-sum estimate `∑_{k≤N} ‖k^{-z}‖ ≤ N` by the sharp
`∑_{1<k≤N} k^{-a} ≤ (N^{1-a}-1)/(1-a)`, uniformly over `sphere c R` with `a := c.re - R ∈ (0,1)`:

  `‖ζ z‖ ≤ (N^{1-a}-1)/(1-a) + 1 + N^{1-a}/(‖c-1‖-R) + (‖c‖+R)·N^{-a}/a`.

With `N ≍ |t|` every term is `≪ |t|^{1-a}`, i.e. log M ≈ (1-a)·log|t| instead of log|t| —
the input for cutting the leading constant of the zero-free region (wiring = next lane-1 brick).
Standard (Titchmarsh §3.5-type); no novelty claimed.  The sphere may straddle Re = 1: the bound
only uses `z.re ≥ a`, monotonicity in the exponent, and `a ≠ 1`.
-/

open Complex Metric Finset MeasureTheory Set

theorem sum_rpow_neg_mono_exp {N : ℕ} {a x : ℝ} (hax : a ≤ x) :
    ∑ k ∈ Ioc 1 N, (k : ℝ) ^ (-x) ≤ ∑ k ∈ Ioc 1 N, (k : ℝ) ^ (-a) := by
  refine Finset.sum_le_sum (fun k hk => ?_)
  have hk1 : 1 < k := (Finset.mem_Ioc.mp hk).1
  have hk1R : (1 : ℝ) ≤ (k : ℝ) := by exact_mod_cast le_of_lt hk1
  exact Real.rpow_le_rpow_of_exponent_le hk1R (by linarith)

theorem norm_riemannZeta_le_sphere_sharpN {c : ℂ} {R : ℝ} {N : ℕ} (hN : 1 ≤ N) (hR : 0 < R)
    (hcR : 0 < c.re - R) (ha1 : c.re - R < 1) (hc1 : R < ‖c - 1‖) {z : ℂ}
    (hz : z ∈ sphere c R) :
    ‖riemannZeta z‖ ≤ ((N : ℝ) ^ (1 - (c.re - R)) - 1) / (1 - (c.re - R)) + 1
      + (N : ℝ) ^ (1 - (c.re - R)) / (‖c - 1‖ - R)
      + (‖c‖ + R) * ((N : ℝ) ^ (-(c.re - R)) / (c.re - R)) := by
  set a : ℝ := c.re - R with ha
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  have hNpos : (0 : ℝ) < (N : ℝ) := by linarith
  have hzc : ‖z - c‖ = R := by rw [← dist_eq_norm]; exact (mem_sphere.mp hz)
  have hzre : a ≤ z.re := by
    have h : |z.re - c.re| ≤ R := by
      rw [← hzc]; simpa [Complex.sub_re] using Complex.abs_re_le_norm (z - c)
    rw [abs_le] at h; rw [ha]; linarith [h.1]
  have hz1 : ‖c - 1‖ - R ≤ ‖z - 1‖ := by
    have htri : ‖c - 1‖ ≤ ‖c - z‖ + ‖z - 1‖ := by
      have : c - 1 = (c - z) + (z - 1) := by ring
      rw [this]; exact norm_add_le _ _
    have hcz : ‖c - z‖ = R := by rw [norm_sub_rev]; exact hzc
    linarith
  have hzn : ‖z‖ ≤ ‖c‖ + R := by
    have : z = c + (z - c) := by ring
    rw [this]; exact le_trans (norm_add_le _ _) (by rw [hzc])
  have hzre0 : 0 < z.re := lt_of_lt_of_le hcR hzre
  have hz1pos : 0 < ‖c - 1‖ - R := by linarith
  have hzne1 : z ≠ 1 := by
    intro h; rw [h] at hz1; simp at hz1; linarith
  set A : ℂ := ∑ k ∈ Finset.Ioc 1 N, (k : ℂ) ^ (-z) with hA
  set B : ℂ := (N : ℂ) ^ (1 - z) / (z - 1) with hB
  set T : ℂ := ∫ t in Ioi (N : ℝ), ((Int.fract t : ℝ) : ℂ) * (t : ℂ) ^ (-(z + 1)) with hT
  rw [riemannZeta_eq_partial_sum_add_tail hN hzre0 hzne1, ← hA, ← hB, ← hT]
  have htri : ‖A + 1 + B - z * T‖ ≤ ‖A‖ + 1 + ‖B‖ + ‖z‖ * ‖T‖ := by
    have h1 : ‖A + 1 + B - z * T‖ ≤ ‖A + 1 + B‖ + ‖z * T‖ := norm_sub_le _ _
    have h2 : ‖A + 1 + B‖ ≤ ‖A + 1‖ + ‖B‖ := norm_add_le _ _
    have h3 : ‖A + 1‖ ≤ ‖A‖ + ‖(1 : ℂ)‖ := norm_add_le _ _
    rw [norm_mul] at h1
    rw [norm_one] at h3
    linarith
  have hAbound : ‖A‖ ≤ ((N : ℝ) ^ (1 - a) - 1) / (1 - a) := by
    refine le_trans norm_sum_cpow_neg_le ?_
    refine le_trans (sum_rpow_neg_mono_exp hzre) ?_
    exact sum_rpow_neg_le hN (le_of_lt hcR) (ne_of_lt ha1)
  have hBnorm : ‖B‖ = (N : ℝ) ^ (1 - z.re) / ‖z - 1‖ := by
    rw [hB, norm_div, ← Complex.ofReal_natCast N, Complex.norm_cpow_eq_rpow_re_of_pos hNpos,
      Complex.sub_re, Complex.one_re]
  have hBbound : ‖B‖ ≤ (N : ℝ) ^ (1 - a) / (‖c - 1‖ - R) := by
    rw [hBnorm]
    have hnum : (N : ℝ) ^ (1 - z.re) ≤ (N : ℝ) ^ (1 - a) :=
      Real.rpow_le_rpow_of_exponent_le hNR (by linarith)
    have hnn : (0 : ℝ) ≤ (N : ℝ) ^ (1 - a) := Real.rpow_nonneg (le_of_lt hNpos) _
    calc (N : ℝ) ^ (1 - z.re) / ‖z - 1‖
        ≤ (N : ℝ) ^ (1 - a) / ‖z - 1‖ :=
          div_le_div_of_nonneg_right hnum (norm_nonneg _)
      _ ≤ (N : ℝ) ^ (1 - a) / (‖c - 1‖ - R) := div_le_div_of_nonneg_left hnn hz1pos hz1
  have hTbound0 : ‖T‖ ≤ (N : ℝ) ^ (-z.re) / z.re :=
    norm_fract_integral_tail_le hzre0 hNR
  have hTbound : ‖T‖ ≤ (N : ℝ) ^ (-a) / a := by
    refine le_trans hTbound0 ?_
    have hnum : (N : ℝ) ^ (-z.re) ≤ (N : ℝ) ^ (-a) :=
      Real.rpow_le_rpow_of_exponent_le hNR (by linarith)
    have hnn : (0 : ℝ) ≤ (N : ℝ) ^ (-a) := Real.rpow_nonneg (le_of_lt hNpos) _
    calc (N : ℝ) ^ (-z.re) / z.re
        ≤ (N : ℝ) ^ (-a) / z.re := div_le_div_of_nonneg_right hnum (le_of_lt hzre0)
      _ ≤ (N : ℝ) ^ (-a) / a := div_le_div_of_nonneg_left hnn hcR hzre
  have hTnn : 0 ≤ ‖T‖ := norm_nonneg _
  have hzT : ‖z‖ * ‖T‖ ≤ (‖c‖ + R) * ((N : ℝ) ^ (-a) / a) :=
    mul_le_mul hzn hTbound hTnn (by linarith [norm_nonneg c])
  calc ‖A + 1 + B - z * T‖
      ≤ ‖A‖ + 1 + ‖B‖ + ‖z‖ * ‖T‖ := htri
    _ ≤ ((N : ℝ) ^ (1 - a) - 1) / (1 - a) + 1 + (N : ℝ) ^ (1 - a) / (‖c - 1‖ - R)
          + (‖c‖ + R) * ((N : ℝ) ^ (-a) / a) := by linarith

/-- Clean form dropping the `-1`: every term is a multiple of `N^{1-a}` or `N^{-a}`. -/
theorem norm_riemannZeta_le_sphere_sharpN_clean {c : ℂ} {R : ℝ} {N : ℕ} (hN : 1 ≤ N)
    (hR : 0 < R) (hcR : 0 < c.re - R) (ha1 : c.re - R < 1) (hc1 : R < ‖c - 1‖) {z : ℂ}
    (hz : z ∈ sphere c R) :
    ‖riemannZeta z‖ ≤ (N : ℝ) ^ (1 - (c.re - R)) * (1 / (1 - (c.re - R)) + 1 / (‖c - 1‖ - R))
      + 1 + (‖c‖ + R) * ((N : ℝ) ^ (-(c.re - R)) / (c.re - R)) := by
  have h := norm_riemannZeta_le_sphere_sharpN hN hR hcR ha1 hc1 hz
  have hpos : 0 < 1 - (c.re - R) := by linarith
  have hdrop : ((N : ℝ) ^ (1 - (c.re - R)) - 1) / (1 - (c.re - R))
      ≤ (N : ℝ) ^ (1 - (c.re - R)) / (1 - (c.re - R)) :=
    div_le_div_of_nonneg_right (by linarith) (le_of_lt hpos)
  have hsplit : (N : ℝ) ^ (1 - (c.re - R)) * (1 / (1 - (c.re - R)) + 1 / (‖c - 1‖ - R))
      = (N : ℝ) ^ (1 - (c.re - R)) / (1 - (c.re - R))
        + (N : ℝ) ^ (1 - (c.re - R)) / (‖c - 1‖ - R) := by ring
  rw [hsplit]; linarith

#print axioms norm_riemannZeta_le_sphere_sharpN
#print axioms norm_riemannZeta_le_sphere_sharpN_clean
