import MulLogLeSelf

/-!
# S2949 — tighter closing estimate `32·log((d·x+10)(3+6a+640·log x)) ≤ 40·log x` (d ≤ 15, log x ≥ 256)

Same decomposition as S2930/§1461 (`argument ≤ 16(3+6a+640)·x·log x`), but the two absorbed terms are charged only
4·log x each: `32·log(16(3+6a+640)) ≤ 4·log x` (threshold `8·log(16(3+6a+640)) ≤ log x`) and
`32·log log x ≤ 4·log x` (`mul_log_le_self` with C = 8, needs log x ≥ 256).  Total 4L + 32L + 4L = 40L,
versus 64L — feeding B = C₂ = C₃ = 40·log|γ| into the 5-8-4-1 region.
-/

theorem log_sharp_close_40 {a x d : ℝ} (ha : 0 ≤ a) (hd : 0 ≤ d) (hd15 : d ≤ 15)
    (hxe : Real.exp 256 ≤ x)
    (hthr : 8 * Real.log (16 * (3 + 6 * a + 640)) ≤ Real.log x) :
    32 * Real.log ((d * x + 10) * (3 + 6 * a + 640 * Real.log x)) ≤ 40 * Real.log x := by
  have hx10 : (10 : ℝ) ≤ x := le_trans (by nlinarith [Real.add_one_le_exp (256 : ℝ)]) hxe
  have hx1 : (1 : ℝ) ≤ x := by linarith
  have hxpos : 0 < x := by linarith
  have hL : 256 ≤ Real.log x := by
    calc (256 : ℝ) = Real.log (Real.exp 256) := (Real.log_exp 256).symm
      _ ≤ Real.log x := Real.log_le_log (Real.exp_pos 256) hxe
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
  have hloglog : Real.log (Real.log x) ≤ Real.log x / 8 := by
    have h := mul_log_le_self (C := 8) (x := Real.log x) (by norm_num) (by norm_num; linarith)
    linarith
  linarith [hlog1, hlog2, hthr, hloglog, hL]

#print axioms log_sharp_close_40
