import Mathlib.NumberTheory.Chebyshev
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

set_option maxHeartbeats 1000000

open Real intervalIntegral

/-!
# §1706a — THE Li PARTS IDENTITY: `∫₂ˣ dt/log t − ∫₂ˣ dt/log²t = x/log x − 2/log 2`

Integration by parts for the logarithmic integral: `d/dt (t/log t) = 1/log t − 1/log²t` on `[2,x]`,
so the offset `Li(x) = ∫₂ˣ dt/log t` differs from `∫₂ˣ dt/log²t` by the boundary term exactly.
This is the bridge between Mathlib's Abel form of `π` (`θ/log + ∫θ/(t log²)`) and `Li(x)`:
substituting `θ = t + err` produces `Li(x) + 2/log 2 + error terms`.
-/

theorem continuousOn_one_div_log {x : ℝ} (hx : 2 ≤ x) :
    ContinuousOn (fun t : ℝ => 1 / Real.log t) (Set.uIcc 2 x) := by
  rw [Set.uIcc_of_le hx]
  intro t ht
  have ht2 : (2 : ℝ) ≤ t := ht.1
  have hlt : 0 < Real.log t := Real.log_pos (by linarith)
  exact (ContinuousAt.div continuousAt_const
    (Real.continuousAt_log (by linarith)) hlt.ne').continuousWithinAt

theorem continuousOn_one_div_log_sq {x : ℝ} (hx : 2 ≤ x) :
    ContinuousOn (fun t : ℝ => 1 / Real.log t ^ 2) (Set.uIcc 2 x) := by
  rw [Set.uIcc_of_le hx]
  intro t ht
  have ht2 : (2 : ℝ) ≤ t := ht.1
  have hlt : 0 < Real.log t := Real.log_pos (by linarith)
  exact (ContinuousAt.div continuousAt_const
    ((Real.continuousAt_log (by linarith)).pow 2) (pow_ne_zero 2 hlt.ne')).continuousWithinAt

theorem li_parts_identity {x : ℝ} (hx : 2 ≤ x) :
    (∫ t in (2:ℝ)..x, 1 / Real.log t) - ∫ t in (2:ℝ)..x, 1 / Real.log t ^ 2
      = x / Real.log x - 2 / Real.log 2 := by
  have hderiv : ∀ t ∈ Set.uIcc (2:ℝ) x,
      HasDerivAt (fun s => s / Real.log s) (1 / Real.log t - 1 / Real.log t ^ 2) t := by
    intro t ht
    rw [Set.uIcc_of_le hx] at ht
    have ht2 : (2 : ℝ) ≤ t := ht.1
    have ht0 : (0 : ℝ) < t := by linarith
    have hlt : 0 < Real.log t := Real.log_pos (by linarith)
    have hd : HasDerivAt (fun s => s / Real.log s)
        ((1 * Real.log t - t * t⁻¹) / (Real.log t) ^ 2) t :=
      HasDerivAt.div (hasDerivAt_id t) (Real.hasDerivAt_log ht0.ne') hlt.ne'
    have heq : (1 * Real.log t - t * t⁻¹) / (Real.log t) ^ 2
        = 1 / Real.log t - 1 / Real.log t ^ 2 := by
      rw [mul_inv_cancel₀ ht0.ne']
      field_simp
    rwa [heq] at hd
  have hint1 : IntervalIntegrable (fun t : ℝ => 1 / Real.log t) MeasureTheory.volume 2 x :=
    (continuousOn_one_div_log hx).intervalIntegrable
  have hint2 : IntervalIntegrable (fun t : ℝ => 1 / Real.log t ^ 2) MeasureTheory.volume 2 x :=
    (continuousOn_one_div_log_sq hx).intervalIntegrable
  have hftc := intervalIntegral.integral_eq_sub_of_hasDerivAt hderiv (hint1.sub hint2)
  rw [← intervalIntegral.integral_sub hint1 hint2]
  exact hftc

#print axioms li_parts_identity
