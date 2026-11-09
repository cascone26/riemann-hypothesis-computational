import GaussWeight

/-!
# S2988 — Gaussian moment bound WITH the log factor (lane 1: three-lines edges, exact exponent + log)

The S2982 edges carry log t: |E(it)| ≪ (1+|t|)^{3/2}(log(1+|t|) + c), |E(1+it)| ≪ (1+|t|)(log(1+|t|) + c).
* `log_one_add_abs_le` — log(1+|t|) ≤ log(1+|T|) + |t − T|.
* `gauss_poly_log_weight` — p, c ≥ 0 ⇒
  (1+|t|)^p (log(1+|t|) + c) e^{−(t−T)²} ≤ (1+|T|)^p (log(1+|T|) + c + 1) e^{(p+1)²/4}   (all real t, T).
Same exponent, same log — the localised edge constants a, b of S2984 are now explicit in T.
-/

namespace GaussWeightLog

theorem log_one_add_abs_le (t T : ℝ) : Real.log (1 + |t|) ≤ Real.log (1 + |T|) + |t - T| := by
  have h1 := GaussWeight.one_add_abs_le_mul t T
  have hp1 : 0 < 1 + |t| := by positivity
  have hpT : 0 < 1 + |T| := by positivity
  have hpu : 0 < 1 + |t - T| := by positivity
  have h2 : Real.log (1 + |t|) ≤ Real.log ((1 + |T|) * (1 + |t - T|)) := Real.log_le_log hp1 h1
  rw [Real.log_mul hpT.ne' hpu.ne'] at h2
  have h3 : Real.log (1 + |t - T|) ≤ |t - T| := by
    have := Real.log_le_sub_one_of_pos hpu; linarith
  linarith

theorem gauss_poly_log_weight (t T : ℝ) {p c : ℝ} (hp : 0 ≤ p) (hc : 0 ≤ c) :
    (1 + |t|) ^ p * (Real.log (1 + |t|) + c) * Real.exp (-(t - T) ^ 2)
      ≤ (1 + |T|) ^ p * (Real.log (1 + |T|) + c + 1) * Real.exp ((p + 1) ^ 2 / 4) := by
  set u := |t - T| with hu
  have hu0 : 0 ≤ u := abs_nonneg _
  have hsq : (t - T) ^ 2 = u ^ 2 := by rw [hu, sq_abs]
  set LT := Real.log (1 + |T|) with hLT
  have hLT0 : 0 ≤ LT := Real.log_nonneg (by linarith [abs_nonneg T])
  have hLt0 : 0 ≤ Real.log (1 + |t|) := Real.log_nonneg (by linarith [abs_nonneg t])
  have hlog : Real.log (1 + |t|) + c ≤ (LT + c + 1) * (1 + u) := by
    have := log_one_add_abs_le t T
    nlinarith
  have hA : (1 + |t|) ^ p ≤ (1 + |T|) ^ p * (1 + u) ^ p := by
    have := Real.rpow_le_rpow (by positivity) (GaussWeight.one_add_abs_le_mul t T) hp
    rwa [Real.mul_rpow (by positivity) (by positivity)] at this
  have hpow : (1 + u) ^ p * (1 + u) = (1 + u) ^ (p + 1) := by
    rw [Real.rpow_add_one (by positivity : (1 + u) ≠ 0)]
  have hB := GaussWeight.one_add_rpow_le_exp hu0 (show 0 ≤ p + 1 by linarith)
  have hE : Real.exp ((p + 1) * u) * Real.exp (-(t - T) ^ 2) ≤ Real.exp ((p + 1) ^ 2 / 4) := by
    rw [← Real.exp_add, hsq]
    apply Real.exp_le_exp.mpr
    nlinarith [sq_nonneg (u - (p + 1) / 2)]
  have hT0 : 0 ≤ (1 + |T|) ^ p := by positivity
  have hg0 : 0 ≤ Real.exp (-(t - T) ^ 2) := (Real.exp_pos _).le
  have hup0 : 0 ≤ (1 + u) ^ p := by positivity
  have hLc : 0 ≤ LT + c + 1 := by linarith
  have hlogt0 : 0 ≤ Real.log (1 + |t|) + c := by linarith
  calc (1 + |t|) ^ p * (Real.log (1 + |t|) + c) * Real.exp (-(t - T) ^ 2)
      ≤ ((1 + |T|) ^ p * (1 + u) ^ p) * ((LT + c + 1) * (1 + u)) * Real.exp (-(t - T) ^ 2) := by
        gcongr
    _ = (1 + |T|) ^ p * (LT + c + 1) * ((1 + u) ^ p * (1 + u)) * Real.exp (-(t - T) ^ 2) := by ring
    _ = (1 + |T|) ^ p * (LT + c + 1) * ((1 + u) ^ (p + 1) * Real.exp (-(t - T) ^ 2)) := by
        rw [hpow]; ring
    _ ≤ (1 + |T|) ^ p * (LT + c + 1) * (Real.exp ((p + 1) * u) * Real.exp (-(t - T) ^ 2)) := by
        gcongr
    _ ≤ (1 + |T|) ^ p * (LT + c + 1) * Real.exp ((p + 1) ^ 2 / 4) := by
        gcongr

end GaussWeightLog

#print axioms GaussWeightLog.log_one_add_abs_le
#print axioms GaussWeightLog.gauss_poly_log_weight
