import MulLogLeSelf

set_option maxHeartbeats 1000000

/-!
# §1565 — sharp slack closing estimate `32·log((d·x+10)(3+6a+640·log x)) ≤ 64·log x - k`

Slack variant of ζ §1461 (`log_sharp_close`), for the QUADRATIC sharp capstone: the χ² correction
`k = ∑_{p|N} log p` is absorbed by tightening the threshold to `2·log(16(3+6a+640)) + k/16 ≤ log x`, giving
`32·log(base) ≤ 16·log x - k`; with `32·log log x ≤ 16·log x` (§1445, C=2) the total is `≤ 64·log x - k`.
EXIT0 / std-axioms.
-/

theorem log_sharp_close_slack {a x d k : ℝ} (ha : 0 ≤ a) (hd : 0 ≤ d) (hd6 : d ≤ 6)
    (hxe : Real.exp 16 ≤ x)
    (hthr : 2 * Real.log (16 * (3 + 6 * a + 640)) + k / 16 ≤ Real.log x) :
    32 * Real.log ((d * x + 10) * (3 + 6 * a + 640 * Real.log x)) ≤ 64 * Real.log x - k := by
  have hx1 : (1 : ℝ) ≤ x := le_trans (by nlinarith [Real.add_one_le_exp (16 : ℝ)]) hxe
  have hxpos : 0 < x := by linarith
  have hL : 16 ≤ Real.log x := by
    calc (16 : ℝ) = Real.log (Real.exp 16) := (Real.log_exp 16).symm
      _ ≤ Real.log x := Real.log_le_log (Real.exp_pos 16) hxe
  have hLpos : 0 < Real.log x := by linarith
  have hbase : 0 < 16 * (3 + 6 * a + 640) := by nlinarith [ha]
  have hprod : (d * x + 10) * (3 + 6 * a + 640 * Real.log x)
      ≤ 16 * (3 + 6 * a + 640) * (x * Real.log x) := by
    nlinarith [ha, hd, hd6, hx1, hL, hLpos, hxpos,
      mul_nonneg (show (0:ℝ) ≤ (16 - d) * x - 10 by nlinarith [hx1]) (show (0:ℝ) ≤ Real.log x by linarith),
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

#print axioms log_sharp_close_slack
