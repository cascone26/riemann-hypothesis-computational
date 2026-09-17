import MulLogLeSelf

set_option maxHeartbeats 1000000

/-!
# §1548 — L closing estimate `32·log(N'·(e·x+6)·(3+6a+10x)) ≤ x` for `e ≤ 16`

L-analogue of ζ §1459, for the L clean bound's first factor `N'·(e·x+6)` (modulus `N' ≥ 1`, coefficient `e` for
both `e=2` at iγ and `e=4` at 2iγ; no `+10` pole constant).  Bound the argument by `(N'(e+6)(16+6a))·x²` (via
`e·x+6 ≤ (e+6)x`, `3+6a+10x ≤ (16+6a)x` for `x ≥ 1`), so `32·log(…) ≤ 32·log(N'(e+6)(16+6a)) + 64·log x
≤ x/4 + x/2 ≤ x` (§1445 at `x ≥ 65536`, base threshold `hx2`).  Closes both `hBclean` (e=2) and `hCclean` (e=4)
for the choice `B=C=|γ|`.  EXIT0 / std-axioms.
-/

theorem log32_close_LFunction {N' a x e : ℝ} (hN' : 1 ≤ N') (ha : 0 ≤ a) (he0 : 0 ≤ e)
    (hx1 : 65536 ≤ x) (hx2 : 128 * Real.log (N' * (e + 6) * (16 + 6 * a)) ≤ x) :
    32 * Real.log (N' * (e * x + 6) * (3 + 6 * a + 10 * x)) ≤ x := by
  have hxpos : 0 < x := by linarith
  have hx1' : (1 : ℝ) ≤ x := by linarith
  have hN'0 : 0 < N' := by linarith
  have he6 : 0 < e + 6 := by linarith
  have hbase : 0 < N' * (e + 6) * (16 + 6 * a) := by positivity
  -- inner product bound (no N'): (e·x+6)(3+6a+10x) ≤ (e+6)(16+6a)x²
  have hinner : (e * x + 6) * (3 + 6 * a + 10 * x) ≤ (e + 6) * (16 + 6 * a) * x ^ 2 := by
    nlinarith [ha, hx1', hxpos, he0,
      mul_nonneg ha (mul_nonneg hxpos.le hxpos.le),
      mul_nonneg (show (0:ℝ) ≤ x - 1 by linarith) hxpos.le,
      mul_nonneg ha (show (0:ℝ) ≤ x - 1 by linarith),
      mul_nonneg he0 (mul_nonneg (show (0:ℝ) ≤ x - 1 by linarith) hxpos.le),
      mul_nonneg (mul_nonneg he0 ha) (mul_nonneg hxpos.le hxpos.le)]
  -- full product bound with N'
  have hprod : N' * (e * x + 6) * (3 + 6 * a + 10 * x)
      ≤ (N' * (e + 6) * (16 + 6 * a)) * x ^ 2 := by
    have h1 : N' * ((e * x + 6) * (3 + 6 * a + 10 * x)) ≤ N' * ((e + 6) * (16 + 6 * a) * x ^ 2) :=
      mul_le_mul_of_nonneg_left hinner hN'0.le
    nlinarith [h1]
  have hargpos : 0 < N' * (e * x + 6) * (3 + 6 * a + 10 * x) := by
    have h1 : 0 < e * x + 6 := by positivity
    have h2 : 0 < 3 + 6 * a + 10 * x := by positivity
    positivity
  have hlog : Real.log (N' * (e * x + 6) * (3 + 6 * a + 10 * x))
      ≤ Real.log ((N' * (e + 6) * (16 + 6 * a)) * x ^ 2) := Real.log_le_log hargpos hprod
  have hsplit : Real.log ((N' * (e + 6) * (16 + 6 * a)) * x ^ 2)
      = Real.log (N' * (e + 6) * (16 + 6 * a)) + 2 * Real.log x := by
    rw [Real.log_mul (ne_of_gt hbase) (by positivity), Real.log_pow]; push_cast; ring
  have h445 : 128 * Real.log x ≤ x := mul_log_le_self (by norm_num) (by norm_num; linarith)
  calc 32 * Real.log (N' * (e * x + 6) * (3 + 6 * a + 10 * x))
      ≤ 32 * Real.log ((N' * (e + 6) * (16 + 6 * a)) * x ^ 2) := by linarith [hlog]
    _ = 32 * Real.log (N' * (e + 6) * (16 + 6 * a)) + 64 * Real.log x := by rw [hsplit]; ring
    _ ≤ x := by linarith [hx2, h445]

#print axioms log32_close_LFunction
