import MulLogLeSelf

set_option maxHeartbeats 1000000

/-!
# S2930 — wide sharp closing estimate `32·log((d·x+10)(3+6a+640·log x)) ≤ 64·log x` for `d ≤ 15`

Widens ζ §1461 (`log_sharp_close`, `d ≤ 6`, covering `iγ` (d=3) and `2iγ` (d=6)) to `d ∈ [0,15]`, so the
same log absorption serves the higher harmonics `3iγ` (d=9), `4iγ` (d=12), `5iγ` (d=15) needed by
higher-degree nonnegative cosine polynomials (5-8-4-1 uses up to `3iγ`).  Only change: the
`(16-d)·x - 10 ≥ 0` slack now uses `x ≥ e^16 ≥ 10` instead of `x ≥ 1`.
-/

theorem log_sharp_close_wide {a x d : ℝ} (ha : 0 ≤ a) (hd : 0 ≤ d) (hd15 : d ≤ 15)
    (hxe : Real.exp 16 ≤ x)
    (hthr : 2 * Real.log (16 * (3 + 6 * a + 640)) ≤ Real.log x) :
    32 * Real.log ((d * x + 10) * (3 + 6 * a + 640 * Real.log x)) ≤ 64 * Real.log x := by
  have hx10 : (10 : ℝ) ≤ x := le_trans (by nlinarith [Real.add_one_le_exp (16 : ℝ)]) hxe
  have hx1 : (1 : ℝ) ≤ x := by linarith
  have hxpos : 0 < x := by linarith
  have hL : 16 ≤ Real.log x := by
    calc (16 : ℝ) = Real.log (Real.exp 16) := (Real.log_exp 16).symm
      _ ≤ Real.log x := Real.log_le_log (Real.exp_pos 16) hxe
  have hLpos : 0 < Real.log x := by linarith
  have hbase : 0 < 16 * (3 + 6 * a + 640) := by nlinarith [ha]
  have hslack : (0:ℝ) ≤ (16 - d) * x - 10 := by nlinarith [hx10, hd15]
  have hprod : (d * x + 10) * (3 + 6 * a + 640 * Real.log x)
      ≤ 16 * (3 + 6 * a + 640) * (x * Real.log x) := by
    nlinarith [ha, hd, hd15, hx1, hL, hLpos, hxpos,
      mul_nonneg hslack (show (0:ℝ) ≤ Real.log x by linarith),
      mul_nonneg hslack (show (0:ℝ) ≤ 3 + 6 * a by linarith),
      mul_nonneg (mul_nonneg ha (show (0:ℝ) ≤ x by linarith)) (show (0:ℝ) ≤ Real.log x - 1 by linarith),
      mul_nonneg (show (0:ℝ) ≤ x by linarith) (show (0:ℝ) ≤ Real.log x - 1 by linarith),
      mul_nonneg (show (0:ℝ) ≤ Real.log x by linarith) (show (0:ℝ) ≤ x by linarith)]
  have hargpos : 0 < (d * x + 10) * (3 + 6 * a + 640 * Real.log x) :=
    mul_pos (by nlinarith [mul_nonneg hd (le_of_lt hxpos)]) (by nlinarith [ha, hLpos])
  have hlog1 : Real.log ((d * x + 10) * (3 + 6 * a + 640 * Real.log x))
      ≤ Real.log (16 * (3 + 6 * a + 640) * (x * Real.log x)) := Real.log_le_log hargpos hprod
  have hlog2 : Real.log (16 * (3 + 6 * a + 640) * (x * Real.log x))
      = Real.log (16 * (3 + 6 * a + 640)) + Real.log x + Real.log (Real.log x) := by
    rw [Real.log_mul (ne_of_gt hbase) (by positivity), Real.log_mul (ne_of_gt hxpos) (ne_of_gt hLpos)]
    ring
  have hloglog : Real.log (Real.log x) ≤ Real.log x / 2 := by
    have h := mul_log_le_self (C := 2) (x := Real.log x) (by norm_num) (by norm_num; linarith)
    linarith
  linarith [hlog1, hlog2, hthr, hloglog, hL]

/-- The `3iγ` instance (d = 9) used by the 5-8-4-1 polynomial. -/
theorem log_sharp_close_3 {a x : ℝ} (ha : 0 ≤ a) (hxe : Real.exp 16 ≤ x)
    (hthr : 2 * Real.log (16 * (3 + 6 * a + 640)) ≤ Real.log x) :
    32 * Real.log ((9 * x + 10) * (3 + 6 * a + 640 * Real.log x)) ≤ 64 * Real.log x :=
  log_sharp_close_wide ha (by norm_num) (by norm_num) hxe hthr

#print axioms log_sharp_close_wide
#print axioms log_sharp_close_3
