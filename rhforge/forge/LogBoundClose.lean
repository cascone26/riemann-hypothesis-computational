import MulLogLeSelf

/-!
# §1458 — the closing transcendental estimate `32·log((3x+10)(3+6a+10x)) ≤ x`

For `a ≥ 0`, `x ≥ 65536`, and `128·log(13+6a) ≤ x`: bound the argument by `((13+6a)x)²`, so
`32·log(…) ≤ 64·log((13+6a)x) = 64log(13+6a) + 64log x ≤ x/2 + x/2 = x` (the last `64log x ≤ x/2` from §1445 with
`C=128` at `x ≥ (256)² = 65536`).  This closes `hBclean` for the choice `B = |γ|`, `K = 3A₀+5|γ|`
(`3+2K = 3+6A₀+10|γ|`), completing the constant-choice fixed point.  EXIT0 / std-axioms.
-/

theorem log_bound_close {a x : ℝ} (ha : 0 ≤ a) (hx1 : 65536 ≤ x) (hx2 : 128 * Real.log (13 + 6 * a) ≤ x) :
    32 * Real.log ((3 * x + 10) * (3 + 6 * a + 10 * x)) ≤ x := by
  have hxpos : 0 < x := by linarith
  have hbase : 0 < 13 + 6 * a := by linarith
  have hbx : 0 < (13 + 6 * a) * x := by positivity
  have hprod : (3 * x + 10) * (3 + 6 * a + 10 * x) ≤ ((13 + 6 * a) * x) ^ 2 := by
    nlinarith [ha, hx1, mul_nonneg (show (0:ℝ) ≤ x - 65536 by linarith) (show (0:ℝ) ≤ x by linarith),
      mul_nonneg (mul_nonneg ha (show (0:ℝ) ≤ x - 65536 by linarith)) (show (0:ℝ) ≤ x by linarith),
      mul_nonneg (mul_nonneg ha ha) (mul_nonneg (le_of_lt hxpos) (le_of_lt hxpos)),
      mul_nonneg ha (le_of_lt hxpos)]
  have hargpos : 0 < (3 * x + 10) * (3 + 6 * a + 10 * x) := by positivity
  have hpoweq : Real.log (((13 + 6 * a) * x) ^ 2) = 2 * Real.log ((13 + 6 * a) * x) := by
    rw [Real.log_pow]; norm_num
  have hlog : Real.log ((3 * x + 10) * (3 + 6 * a + 10 * x)) ≤ 2 * Real.log ((13 + 6 * a) * x) := by
    rw [← hpoweq]; exact Real.log_le_log hargpos hprod
  have hsplit : Real.log ((13 + 6 * a) * x) = Real.log (13 + 6 * a) + Real.log x :=
    Real.log_mul (ne_of_gt hbase) (ne_of_gt hxpos)
  have h445 : 128 * Real.log x ≤ x := mul_log_le_self (by norm_num) (by norm_num; linarith)
  calc 32 * Real.log ((3 * x + 10) * (3 + 6 * a + 10 * x))
      ≤ 32 * (2 * Real.log ((13 + 6 * a) * x)) := by linarith [hlog]
    _ = 64 * Real.log (13 + 6 * a) + 64 * Real.log x := by rw [hsplit]; ring
    _ ≤ x := by linarith [hx2, h445]

#print axioms log_bound_close
