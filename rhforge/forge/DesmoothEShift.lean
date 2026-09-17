import Mathlib.Analysis.SpecialFunctions.Pow.Real

set_option maxHeartbeats 1000000

/-!
# §1704b — E-SHIFT FACTS: `E(x+h) ≤ E(x)` and `E(x-h) ≤ e^{c√log 2}·E(x)` for `h ≤ x/2`

`E(x) = e^{-c√log x}` is decreasing, and shifting DOWN by `h ≤ x/2` costs at most the fixed factor
`e^{c√log 2}`: `log x ≤ log(2(x-h)) = log(x-h) + log 2`, and `√` is subadditive (proved inline —
Mathlib lacks it), so `√log x ≤ √log(x-h) + √log 2`.  These let the de-smoothing squeeze evaluate
all three `Ψ`-errors (at `x-h`, `x`, `x+h`) against the single scale `E(x)`.
-/

theorem sqrt_add_le_sqrt_add_sqrt {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) :
    Real.sqrt (a + b) ≤ Real.sqrt a + Real.sqrt b := by
  have h1 : a + b ≤ (Real.sqrt a + Real.sqrt b) ^ 2 := by
    have hA := Real.sq_sqrt ha
    have hB := Real.sq_sqrt hb
    nlinarith [mul_nonneg (Real.sqrt_nonneg a) (Real.sqrt_nonneg b)]
  have h2 := Real.sqrt_le_sqrt h1
  rwa [Real.sqrt_sq (by positivity)] at h2

theorem desmooth_E_shift {c x h : ℝ} (hc : 0 < c) (hx2 : 2 ≤ x) (hh0 : 0 < h)
    (hhx : h ≤ x / 2) :
    Real.exp (-(c * Real.sqrt (Real.log (x + h))))
        ≤ Real.exp (-(c * Real.sqrt (Real.log x)))
    ∧ Real.exp (-(c * Real.sqrt (Real.log (x - h))))
        ≤ Real.exp (c * Real.sqrt (Real.log 2))
            * Real.exp (-(c * Real.sqrt (Real.log x))) := by
  have hx0 : (0 : ℝ) < x := by linarith
  have hxh1 : (1 : ℝ) ≤ x - h := by linarith
  have hxh0 : (0 : ℝ) < x - h := by linarith
  have hlogxh : (0 : ℝ) ≤ Real.log (x - h) := Real.log_nonneg hxh1
  have hlog2 : (0 : ℝ) ≤ Real.log 2 := Real.log_nonneg (by norm_num)
  constructor
  · apply Real.exp_le_exp.mpr
    have h1 : Real.log x ≤ Real.log (x + h) := Real.log_le_log hx0 (by linarith)
    have h2 := Real.sqrt_le_sqrt h1
    nlinarith [h2]
  · rw [← Real.exp_add]
    apply Real.exp_le_exp.mpr
    -- √log x ≤ √log(x-h) + √log 2
    have hlogle : Real.log x ≤ Real.log (x - h) + Real.log 2 := by
      have h1 : x ≤ (x - h) * 2 := by linarith
      calc Real.log x ≤ Real.log ((x - h) * 2) := Real.log_le_log hx0 h1
        _ = Real.log (x - h) + Real.log 2 := Real.log_mul hxh0.ne' (by norm_num)
    have hsq : Real.sqrt (Real.log x)
        ≤ Real.sqrt (Real.log (x - h)) + Real.sqrt (Real.log 2) := by
      calc Real.sqrt (Real.log x)
          ≤ Real.sqrt (Real.log (x - h) + Real.log 2) := Real.sqrt_le_sqrt hlogle
        _ ≤ Real.sqrt (Real.log (x - h)) + Real.sqrt (Real.log 2) :=
            sqrt_add_le_sqrt_add_sqrt hlogxh hlog2
    nlinarith [hsq]

#print axioms sqrt_add_le_sqrt_add_sqrt
#print axioms desmooth_E_shift
