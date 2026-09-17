import PNTPiAsymptotic

set_option maxHeartbeats 1000000

open Chebyshev Real Filter
open scoped Nat.Prime Topology

/-!
# §1709 ★★★ — THE CLASSIC LIMIT FORM: `π(x)·log x / x → 1`

Capstone corollary: dividing §1708's `|π(x) − x/log x| ≤ C·x/log²x` by `x/log x` gives
`|π(x)·log x/x − 1| ≤ C/log x → 0`.  The prime number theorem in the form Hadamard and
de la Vallée Poussin proved in 1896, as a `Filter.Tendsto`, with the whole riemann9 tower
underneath it.
-/

theorem pnt_ratio_tendsto :
    Tendsto (fun x : ℝ => (π ⌊x⌋₊ : ℝ) * Real.log x / x) atTop (𝓝 1) := by
  obtain ⟨C, x₀, hC, hx₀, hbound⟩ := pnt_pi_asymptotic
  have hzero : Tendsto (fun x : ℝ => (π ⌊x⌋₊ : ℝ) * Real.log x / x - 1) atTop (𝓝 0) := by
    apply squeeze_zero_norm' (a := fun x : ℝ => C / Real.log x)
    · -- eventual pointwise bound
      filter_upwards [eventually_ge_atTop (max x₀ (Real.exp 1))] with x hx
      have hxx₀ : x₀ ≤ x := le_trans (le_max_left _ _) hx
      have hxe : Real.exp 1 ≤ x := le_trans (le_max_right _ _) hx
      have hx0 : (0 : ℝ) < x := lt_of_lt_of_le (Real.exp_pos _) hxe
      have hL1 : (1 : ℝ) ≤ Real.log x := (Real.le_log_iff_exp_le hx0).mpr hxe
      have hL0 : (0 : ℝ) < Real.log x := by linarith
      have hb := hbound x hxx₀
      have heq : (π ⌊x⌋₊ : ℝ) * Real.log x / x - 1
          = ((π ⌊x⌋₊ : ℝ) - x / Real.log x) * (Real.log x / x) := by
        field_simp
      rw [Real.norm_eq_abs, heq, abs_mul,
        abs_of_pos (by positivity : (0:ℝ) < Real.log x / x)]
      have hstep : |(π ⌊x⌋₊ : ℝ) - x / Real.log x| * (Real.log x / x)
          ≤ (C * x / Real.log x ^ 2) * (Real.log x / x) :=
        mul_le_mul_of_nonneg_right hb (by positivity)
      have hcollapse : (C * x / Real.log x ^ 2) * (Real.log x / x) = C / Real.log x := by
        field_simp
      rw [hcollapse] at hstep
      exact hstep
    · -- C / log x → 0
      have h1 : Tendsto Real.log atTop atTop := Real.tendsto_log_atTop
      have h2 : Tendsto (fun x : ℝ => (Real.log x)⁻¹) atTop (𝓝 0) :=
        h1.inv_tendsto_atTop
      have h3 := h2.const_mul C
      simp only [mul_zero] at h3
      convert h3 using 2 with x
      rw [div_eq_mul_inv]
  have hfinal := hzero.add_const 1
  simp only [zero_add] at hfinal
  convert hfinal using 2 with x
  ring

#print axioms pnt_ratio_tendsto
