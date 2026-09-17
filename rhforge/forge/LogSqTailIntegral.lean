import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Tactic

set_option maxHeartbeats 2000000

open intervalIntegral Real

/-!
# §1693 — the CONVERGENT TAIL: `∫_a^b log²(t+2)/t² dt ≤ 96/√a`

Why the left-line error integral stays BOUNDED uniformly in `T`: the §1692 large-`t` integrand
`log²(t+2)/t²` is dominated by `28·t^{-3/2}` (via `log y ≤ 4·y^{1/4}`, from
`log ≤ x-1` applied to the fourth root), whose tail integral is `≤ 56/√a`. Purely real,
absolute constants.
-/

theorem log_le_four_mul_rpow {y : ℝ} (hy : 1 ≤ y) : Real.log y ≤ 4 * y ^ (1/4 : ℝ) := by
  have hypos : (0 : ℝ) < y := by linarith
  have h1 : Real.log y = 4 * Real.log (y ^ (1/4 : ℝ)) := by
    rw [Real.log_rpow hypos]
    ring
  have h2 : Real.log (y ^ (1/4 : ℝ)) ≤ y ^ (1/4 : ℝ) - 1 :=
    Real.log_le_sub_one_of_pos (Real.rpow_pos_of_pos hypos _)
  have h3 : (0 : ℝ) ≤ y ^ (1/4 : ℝ) := (Real.rpow_pos_of_pos hypos _).le
  rw [h1]
  nlinarith

theorem integral_log_sq_div_sq_le {a b : ℝ} (ha : 1 ≤ a) (hab : a ≤ b) :
    (∫ t in a..b, (Real.log (t + 2)) ^ 2 / t ^ 2) ≤ 96 / Real.sqrt a := by
  have hapos : (0 : ℝ) < a := by linarith
  -- pointwise domination by 48·t^(-3/2)
  have hpt : ∀ t ∈ Set.Icc a b, (Real.log (t + 2)) ^ 2 / t ^ 2 ≤ 48 * t ^ (-(3/2) : ℝ) := by
    intro t ht
    have hta : a ≤ t := ht.1
    have ht1 : (1 : ℝ) ≤ t := le_trans ha hta
    have htpos : (0 : ℝ) < t := by linarith
    have hlog : Real.log (t + 2) ≤ 4 * (t + 2) ^ (1/4 : ℝ) :=
      log_le_four_mul_rpow (by linarith)
    have hlognn : 0 ≤ Real.log (t + 2) := Real.log_nonneg (by linarith)
    have hsq : (Real.log (t + 2)) ^ 2 ≤ 16 * (t + 2) ^ (1/2 : ℝ) := by
      have h1 : (Real.log (t + 2)) ^ 2 ≤ (4 * (t + 2) ^ (1/4 : ℝ)) ^ 2 := by
        apply sq_le_sq' (by nlinarith [Real.rpow_pos_of_pos (show (0:ℝ) < t + 2 by linarith) (1/4 : ℝ)]) hlog
      calc (Real.log (t + 2)) ^ 2 ≤ (4 * (t + 2) ^ (1/4 : ℝ)) ^ 2 := h1
        _ = 16 * ((t + 2) ^ (1/4 : ℝ)) ^ 2 := by ring
        _ = 16 * (t + 2) ^ (1/2 : ℝ) := by
            rw [← Real.rpow_natCast ((t + 2) ^ (1/4 : ℝ)) 2, ← Real.rpow_mul (by linarith)]
            norm_num
    have h32 : (t + 2) ^ (1/2 : ℝ) ≤ (3 * t) ^ (1/2 : ℝ) := by
      apply Real.rpow_le_rpow (by linarith) (by linarith) (by norm_num)
    have h3t : (3 * t) ^ (1/2 : ℝ) = Real.sqrt 3 * t ^ (1/2 : ℝ) := by
      rw [Real.mul_rpow (by norm_num) htpos.le, Real.sqrt_eq_rpow]
    have hsqrt3 : Real.sqrt 3 ≤ 2 := by
      nlinarith [Real.sq_sqrt (show (0:ℝ) ≤ 3 by norm_num), Real.sqrt_nonneg 3]
    calc (Real.log (t + 2)) ^ 2 / t ^ 2
        ≤ 16 * (t + 2) ^ (1/2 : ℝ) / t ^ 2 := by gcongr
      _ ≤ 16 * (Real.sqrt 3 * t ^ (1/2 : ℝ)) / t ^ 2 := by
          rw [← h3t]
          gcongr
      _ ≤ 16 * (2 * t ^ (1/2 : ℝ)) / t ^ 2 := by
          gcongr
      _ = 32 * (t ^ (1/2 : ℝ) / t ^ 2) := by ring
      _ = 32 * t ^ (-(3/2) : ℝ) := by
          congr 1
          rw [← Real.rpow_natCast t 2, ← Real.rpow_sub htpos]
          norm_num
      _ ≤ 48 * t ^ (-(3/2) : ℝ) := by
          have := Real.rpow_nonneg htpos.le (-(3/2) : ℝ)
          linarith
  -- integrate the dominating function
  have hint_dom : (∫ t in a..b, 48 * t ^ (-(3/2) : ℝ)) ≤ 96 / Real.sqrt a := by
    have hne : (-(3/2) : ℝ) ≠ -1 := by norm_num
    have hmem : (0 : ℝ) ∉ Set.uIcc a b := by
      rw [Set.uIcc_of_le hab]
      intro h
      have := (Set.mem_Icc.mp h).1
      linarith
    rw [intervalIntegral.integral_const_mul, integral_rpow (Or.inr ⟨hne, hmem⟩)]
    have hb0 : (0 : ℝ) ≤ b ^ (-(3/2) + 1 : ℝ) :=
      Real.rpow_nonneg (by linarith) _
    have ha' : a ^ (-(3/2) + 1 : ℝ) = 1 / Real.sqrt a := by
      rw [show (-(3/2) + 1 : ℝ) = -(1/2) by norm_num, Real.rpow_neg hapos.le,
        Real.sqrt_eq_rpow]
      norm_num
    have : (b ^ (-(3/2) + 1 : ℝ) - a ^ (-(3/2) + 1 : ℝ)) / (-(3/2) + 1)
        = 2 * (a ^ (-(3/2) + 1 : ℝ) - b ^ (-(3/2) + 1 : ℝ)) := by
      rw [show (-(3/2) + 1 : ℝ) = -(1/2) by norm_num]
      ring
    rw [this, ha']
    have hsqrtpos : 0 < Real.sqrt a := Real.sqrt_pos.mpr hapos
    calc 48 * (2 * (1 / Real.sqrt a - b ^ (-(3/2) + 1 : ℝ)))
        ≤ 48 * (2 * (1 / Real.sqrt a)) := by
          gcongr
          linarith
      _ = 96 / Real.sqrt a := by ring
  -- comparison
  have hcont1 : ContinuousOn (fun t : ℝ => (Real.log (t + 2)) ^ 2 / t ^ 2) (Set.uIcc a b) := by
    rw [Set.uIcc_of_le hab]
    apply ContinuousOn.div
    · apply ContinuousOn.pow
      apply ContinuousOn.log
      · fun_prop
      · intro t ht
        have := ht.1
        intro h
        linarith [le_trans ha this]
    · fun_prop
    · intro t ht
      have h1 : (1 : ℝ) ≤ t := le_trans ha ht.1
      have h0 : (0 : ℝ) < t := by linarith
      positivity
  have hcont2 : ContinuousOn (fun t : ℝ => 48 * t ^ (-(3/2) : ℝ)) (Set.uIcc a b) := by
    rw [Set.uIcc_of_le hab]
    intro t ht
    have h1 : (1 : ℝ) ≤ t := le_trans ha ht.1
    apply ContinuousWithinAt.mul continuousWithinAt_const
    apply ContinuousAt.continuousWithinAt
    exact Real.continuousAt_rpow_const t _
      (Or.inl (by intro h; rw [h] at h1; norm_num at h1))
  calc (∫ t in a..b, (Real.log (t + 2)) ^ 2 / t ^ 2)
      ≤ ∫ t in a..b, 48 * t ^ (-(3/2) : ℝ) :=
        intervalIntegral.integral_mono_on hab hcont1.intervalIntegrable
          hcont2.intervalIntegrable hpt
    _ ≤ 96 / Real.sqrt a := hint_dom

#print axioms integral_log_sq_div_sq_le
