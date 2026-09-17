import MulLogLeSelf

set_option maxHeartbeats 1000000

/-!
# §1550 — sharp closing estimate `32·log(N'·(e·x+6)·(3+6a+640·log x)) ≤ 64·log x`

L-analogue of ζ §1461: the `log|γ|`-scaled version for the SHARP region (`B = 64·log|γ|`, `K ≍ log|γ|`, region
`1-c/log|γ|`), for the L first factor `N'·(e·x+6)` (modulus `N' ≥ 1`, coefficient `e`).  Bound the argument by
`N'(e+6)(3+6a+640)·(x·log x)` (via `e·x+6 ≤ (e+6)x`, `3+6a+640log x ≤ (3+6a+640)log x` for `x,log x ≥ 1`), so
`32·log(…) ≤ 32log(N'(e+6)(3+6a+640)) + 32log x + 32log log x`, absorbing `32log(base) ≤ 16log x` (threshold) and
`32log log x ≤ 16log x` (§1445, `C=2`, `log x ≥ 16`).  Sum `≤ 16L+32L+16L = 64L`.  EXIT0 / std-axioms.
-/

theorem log_sharp_close_LFunction {N' a x e : ℝ} (hN' : 1 ≤ N') (ha : 0 ≤ a) (he0 : 0 ≤ e)
    (hxe : Real.exp 16 ≤ x)
    (hthr : 2 * Real.log (N' * (e + 6) * (3 + 6 * a + 640)) ≤ Real.log x) :
    32 * Real.log (N' * (e * x + 6) * (3 + 6 * a + 640 * Real.log x)) ≤ 64 * Real.log x := by
  have hx1 : (1 : ℝ) ≤ x := le_trans (by nlinarith [Real.add_one_le_exp (16 : ℝ)]) hxe
  have hxpos : 0 < x := by linarith
  have hN'0 : 0 < N' := by linarith
  have he6 : 0 < e + 6 := by linarith
  have hL : 16 ≤ Real.log x := by
    calc (16 : ℝ) = Real.log (Real.exp 16) := (Real.log_exp 16).symm
      _ ≤ Real.log x := Real.log_le_log (Real.exp_pos 16) hxe
  have hLpos : 0 < Real.log x := by linarith
  have hbase : 0 < N' * (e + 6) * (3 + 6 * a + 640) := by positivity
  have hinner : (e * x + 6) * (3 + 6 * a + 640 * Real.log x)
      ≤ (e + 6) * (3 + 6 * a + 640) * (x * Real.log x) := by
    nlinarith [ha, he0, hx1, hL, hLpos, hxpos,
      mul_nonneg he0 (show (0:ℝ) ≤ x - 1 by linarith),
      mul_nonneg (show (0:ℝ) ≤ x - 1 by linarith) (show (0:ℝ) ≤ Real.log x by linarith),
      mul_nonneg ha (show (0:ℝ) ≤ Real.log x - 1 by linarith),
      mul_nonneg (mul_nonneg he0 (show (0:ℝ) ≤ x by linarith)) (show (0:ℝ) ≤ Real.log x - 1 by linarith),
      mul_nonneg (show (0:ℝ) ≤ Real.log x by linarith) (show (0:ℝ) ≤ x - 1 by linarith),
      mul_nonneg (mul_nonneg ha (show (0:ℝ) ≤ x - 1 by linarith)) (show (0:ℝ) ≤ Real.log x by linarith)]
  have hprod : N' * (e * x + 6) * (3 + 6 * a + 640 * Real.log x)
      ≤ (N' * (e + 6) * (3 + 6 * a + 640)) * (x * Real.log x) := by
    have h1 : N' * ((e * x + 6) * (3 + 6 * a + 640 * Real.log x))
        ≤ N' * ((e + 6) * (3 + 6 * a + 640) * (x * Real.log x)) :=
      mul_le_mul_of_nonneg_left hinner hN'0.le
    nlinarith [h1]
  have hargpos : 0 < N' * (e * x + 6) * (3 + 6 * a + 640 * Real.log x) := by
    have h1 : 0 < e * x + 6 := by positivity
    have h2 : 0 < 3 + 6 * a + 640 * Real.log x := by positivity
    positivity
  have hlog1 : Real.log (N' * (e * x + 6) * (3 + 6 * a + 640 * Real.log x))
      ≤ Real.log ((N' * (e + 6) * (3 + 6 * a + 640)) * (x * Real.log x)) :=
    Real.log_le_log hargpos hprod
  have hlog2 : Real.log ((N' * (e + 6) * (3 + 6 * a + 640)) * (x * Real.log x))
      = Real.log (N' * (e + 6) * (3 + 6 * a + 640)) + Real.log x + Real.log (Real.log x) := by
    rw [Real.log_mul (ne_of_gt hbase) (by positivity), Real.log_mul (ne_of_gt hxpos) (ne_of_gt hLpos)]
    ring
  have hloglog : Real.log (Real.log x) ≤ Real.log x / 2 := by
    have h := mul_log_le_self (C := 2) (x := Real.log x) (by norm_num) (by norm_num; linarith)
    linarith
  linarith [hlog1, hlog2, hthr, hloglog, hL]

#print axioms log_sharp_close_LFunction
