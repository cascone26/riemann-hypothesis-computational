import ThetaErrorIntegral

set_option maxHeartbeats 4000000

open Chebyshev Real MeasureTheory intervalIntegral
open scoped Nat.Prime

/-!
# §1706 ★★★★★ — THE PRIME NUMBER THEOREM: `π(x) = Li(x) + O(x·e^{-c√log x})`

The end of the chain.  Mathlib's Abel identity `π⌊x⌋₊ = θx/log x + ∫₂ˣ θt/(t log²t)` plus the
parts identity (§1706a) turns the θ-error (§1705) into
`π⌊x⌋₊ − ∫₂ˣ dt/log t = (θx−x)/log x + ∫₂ˣ (θt−t)/(t log²t) + 2/log 2`, and each right-hand term
is `O(x·e^{-(c/2)√log x})` (§1706b for the integral; the constant absorbs since
`x·e^{-(c/2)√L} ≥ e^{L/2} ≥ 3 ≥ 2/log 2`).  This is the prime number theorem with de la Vallée
Poussin's classical error term, EXIT0 from Mathlib's axioms through the entire riemann9 tower:
zero-free region → winding rectangle → Landau ladder → effective Perron → smoothed PNT →
de-smoothing → θ → π.  (Report-tier honesty: formalizes KNOWN mathematics — this is verified
scaffolding for the analytic-number-theory stack, not progress on RH itself.)
-/

theorem pnt_pi_li_error :
    ∃ C c x₀ : ℝ, 0 < C ∧ 0 < c ∧ 2 ≤ x₀ ∧ ∀ x : ℝ, x₀ ≤ x →
      |(π ⌊x⌋₊ : ℝ) - ∫ t in (2:ℝ)..x, 1 / Real.log t|
        ≤ C * x * Real.exp (-(c * Real.sqrt (Real.log x))) := by
  obtain ⟨C, c, x₀θ, hC, hc, hx₀θ1, hθ⟩ := pnt_theta_error
  refine ⟨8 + 4 * C, c / 2,
    max ((max x₀θ 2) ^ 2) (max 4 (Real.exp (max (c ^ 2) 4))),
    by positivity, by positivity,
    le_trans (by nlinarith [le_max_right x₀θ 2] : (2:ℝ) ≤ (max x₀θ 2) ^ 2)
      (le_max_left _ _), ?_⟩
  intro x hx
  have hxsq : (max x₀θ 2) ^ 2 ≤ x := le_trans (le_max_left _ _) hx
  have hx4 : (4 : ℝ) ≤ x := le_trans (le_trans (le_max_left _ _) (le_max_right _ _)) hx
  have hxc : Real.exp (max (c ^ 2) 4) ≤ x :=
    le_trans (le_trans (le_max_right _ _) (le_max_right _ _)) hx
  have hx0 : (0 : ℝ) < x := by linarith
  have hx2 : (2 : ℝ) ≤ x := by linarith
  have hL : max (c ^ 2) 4 ≤ Real.log x := (Real.le_log_iff_exp_le hx0).mpr hxc
  have hLc : c ^ 2 ≤ Real.log x := le_trans (le_max_left _ _) hL
  have hL4 : (4 : ℝ) ≤ Real.log x := le_trans (le_max_right _ _) hL
  have hL0 : (0 : ℝ) ≤ Real.log x := by linarith
  have hsqLc : c ≤ Real.sqrt (Real.log x) := by
    have h2 := Real.sqrt_le_sqrt hLc
    rwa [Real.sqrt_sq hc.le] at h2
  set E := Real.exp (-(c / 2 * Real.sqrt (Real.log x))) with hEdef
  have hE0 : (0 : ℝ) < E := by rw [hEdef]; exact Real.exp_pos _
  -- the three ingredients
  have hid := primeCounting_eq_theta_div_log_add_integral hx2
  have hparts := li_parts_identity hx2
  have hbound' : ∀ t : ℝ, max x₀θ 2 ≤ t →
      |θ t - t| ≤ C * t * Real.exp (-(c * Real.sqrt (Real.log t))) :=
    fun t ht => hθ t (le_trans (le_max_left _ _) ht)
  have hint := theta_error_integral_bound hC hc (le_max_right x₀θ 2) hbound' hxsq hx4
    (le_trans (Real.exp_le_exp.mpr (le_max_left _ _)) hxc)
  -- linearity: ∫ θ/(t log²) = ∫ g + ∫ 1/log²
  have hIθ : IntervalIntegrable (fun t : ℝ => θ t / (t * Real.log t ^ 2))
      MeasureTheory.volume 2 x := by
    rw [intervalIntegrable_iff, Set.uIoc_of_le hx2]
    exact (integrableOn_theta_div_id_mul_log_sq x).mono_set Set.Ioc_subset_Icc_self
  have hI2 : IntervalIntegrable (fun t : ℝ => 1 / Real.log t ^ 2)
      MeasureTheory.volume 2 x := (continuousOn_one_div_log_sq hx2).intervalIntegrable
  have hlin := intervalIntegral.integral_sub hIθ hI2
  -- the decomposition:  π − Li = (θx−x)/log x + ∫g + 2/log 2
  have hlogx1 : (1 : ℝ) ≤ Real.log x := by linarith
  have hdecomp : (π ⌊x⌋₊ : ℝ) - ∫ t in (2:ℝ)..x, 1 / Real.log t
      = (θ x - x) / Real.log x
        + (∫ t in (2:ℝ)..x, (θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2))
        + 2 / Real.log 2 := by
    have h1 : (θ x - x) / Real.log x = θ x / Real.log x - x / Real.log x := sub_div _ _ _
    rw [hid, h1]
    linarith [hlin, hparts]
  rw [hdecomp]
  -- term bounds
  have hT1 : |(θ x - x) / Real.log x| ≤ C * x * E := by
    have h1 : |(θ x - x) / Real.log x| = |θ x - x| / Real.log x := by
      rw [abs_div, abs_of_pos (by linarith : (0:ℝ) < Real.log x)]
    have h2 : |θ x - x| / Real.log x ≤ |θ x - x| :=
      div_le_self (abs_nonneg _) hlogx1
    have h3 := hθ x (by nlinarith [le_max_left x₀θ 2, le_max_right x₀θ 2])
    have h4 : Real.exp (-(c * Real.sqrt (Real.log x))) ≤ E := by
      rw [hEdef]
      apply Real.exp_le_exp.mpr
      nlinarith [Real.sqrt_nonneg (Real.log x), hc.le]
    have h5 : C * x * Real.exp (-(c * Real.sqrt (Real.log x))) ≤ C * x * E :=
      mul_le_mul_of_nonneg_left h4 (by positivity)
    rw [h1]
    linarith
  have hT3 : 2 / Real.log 2 ≤ x * E := by
    have hlog2 : (0.6931471803 : ℝ) < Real.log 2 := Real.log_two_gt_d9
    have h1 : 2 / Real.log 2 ≤ 3 := by
      rw [div_le_iff₀ (by linarith)]
      linarith
    have h2 : (3 : ℝ) ≤ Real.exp 2 := by
      have := Real.add_one_le_exp (2 : ℝ)
      linarith
    have h3 : Real.exp 2 ≤ Real.exp (Real.log x / 2) :=
      Real.exp_le_exp.mpr (by linarith)
    have h4 : Real.exp (Real.log x / 2) ≤ x * E := by
      rw [hEdef]
      have hkey : Real.log x / 2 ≤ Real.log x + -(c / 2 * Real.sqrt (Real.log x)) := by
        have ha := mul_le_mul_of_nonneg_right hsqLc (Real.sqrt_nonneg (Real.log x))
        have hb : Real.sqrt (Real.log x) * Real.sqrt (Real.log x) = Real.log x :=
          Real.mul_self_sqrt hL0
        nlinarith [ha, hb]
      calc Real.exp (Real.log x / 2)
          ≤ Real.exp (Real.log x + -(c / 2 * Real.sqrt (Real.log x))) :=
            Real.exp_le_exp.mpr hkey
        _ = x * Real.exp (-(c / 2 * Real.sqrt (Real.log x))) := by
            rw [Real.exp_add, Real.exp_log hx0]
    linarith
  -- assemble
  calc |(θ x - x) / Real.log x
        + (∫ t in (2:ℝ)..x, (θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2))
        + 2 / Real.log 2|
      ≤ |(θ x - x) / Real.log x
          + (∫ t in (2:ℝ)..x, (θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2))|
        + |2 / Real.log 2| := abs_add_le _ _
    _ ≤ |(θ x - x) / Real.log x|
        + |∫ t in (2:ℝ)..x, (θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2)|
        + |2 / Real.log 2| := by
        have := abs_add_le ((θ x - x) / Real.log x)
          (∫ t in (2:ℝ)..x, (θ t / (t * Real.log t ^ 2) - 1 / Real.log t ^ 2))
        linarith
    _ ≤ C * x * E + (7 + 3 * C) * x * E + x * E := by
        have habs2 : |(2:ℝ) / Real.log 2| = 2 / Real.log 2 := by
          have hlog2 : (0.6931471803 : ℝ) < Real.log 2 := Real.log_two_gt_d9
          exact abs_of_pos (by positivity)
        rw [habs2]
        have := hint
        rw [← hEdef] at this
        linarith [hT1, hT3, this]
    _ = (8 + 4 * C) * x * E := by ring

#print axioms pnt_pi_li_error
