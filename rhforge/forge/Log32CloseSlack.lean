import MulLogLeSelf

set_option maxHeartbeats 1000000

/-!
# §1563 — slack closing estimate `32·log((d·x+10)(3+6a+10x)) ≤ x - k`

Variant of ζ §1459 with a slack `k ≥ 0` (for the finite-product correction `k = ∑_{p|N} log p ≤ log N` of the
QUADRATIC case).  Same argument, stronger threshold `128·log(16+6a) + 2k ≤ x`: from `32·log(…) ≤ 64log(16+6a)
+ 64log x`, absorb `64log(16+6a) ≤ x/2 - k` (threshold) and `64log x ≤ x/2` (§1445).  Sum `≤ x - k`.
EXIT0 / std-axioms.
-/

theorem log32_close_slack {a x d k : ℝ} (ha : 0 ≤ a) (hd : 0 ≤ d) (hd16 : d ≤ 16)
    (hx1 : 65536 ≤ x) (hx2 : 128 * Real.log (16 + 6 * a) + 2 * k ≤ x) :
    32 * Real.log ((d * x + 10) * (3 + 6 * a + 10 * x)) ≤ x - k := by
  have hxpos : 0 < x := by linarith
  have hbase : 0 < 16 + 6 * a := by linarith
  have hprod : (d * x + 10) * (3 + 6 * a + 10 * x) ≤ ((16 + 6 * a) * x) ^ 2 := by
    nlinarith [ha, hd, hd16, hx1,
      mul_nonneg (show (0:ℝ) ≤ x - 65536 by linarith) (show (0:ℝ) ≤ x by linarith),
      mul_nonneg (mul_nonneg ha (show (0:ℝ) ≤ x - 65536 by linarith)) (show (0:ℝ) ≤ x by linarith),
      mul_nonneg (mul_nonneg ha ha) (mul_nonneg (le_of_lt hxpos) (le_of_lt hxpos)),
      mul_nonneg (show (0:ℝ) ≤ 16 - d by linarith) (mul_nonneg (le_of_lt hxpos) (le_of_lt hxpos)),
      mul_nonneg ha (le_of_lt hxpos)]
  have hargpos : 0 < (d * x + 10) * (3 + 6 * a + 10 * x) :=
    mul_pos (by nlinarith [mul_nonneg hd (le_of_lt hxpos)]) (by nlinarith [ha, hxpos])
  have hpoweq : Real.log (((16 + 6 * a) * x) ^ 2) = 2 * Real.log ((16 + 6 * a) * x) := by
    rw [Real.log_pow]; norm_num
  have hlog : Real.log ((d * x + 10) * (3 + 6 * a + 10 * x)) ≤ 2 * Real.log ((16 + 6 * a) * x) := by
    rw [← hpoweq]; exact Real.log_le_log hargpos hprod
  have hsplit : Real.log ((16 + 6 * a) * x) = Real.log (16 + 6 * a) + Real.log x :=
    Real.log_mul (ne_of_gt hbase) (ne_of_gt hxpos)
  have h445 : 128 * Real.log x ≤ x := mul_log_le_self (by norm_num) (by norm_num; linarith)
  calc 32 * Real.log ((d * x + 10) * (3 + 6 * a + 10 * x))
      ≤ 32 * (2 * Real.log ((16 + 6 * a) * x)) := by linarith [hlog]
    _ = 64 * Real.log (16 + 6 * a) + 64 * Real.log x := by rw [hsplit]; ring
    _ ≤ x - k := by linarith [hx2, h445]

#print axioms log32_close_slack
