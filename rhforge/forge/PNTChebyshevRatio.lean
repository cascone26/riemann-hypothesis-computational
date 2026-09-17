import PNTThetaError

set_option maxHeartbeats 1000000

open Chebyshev Real Filter
open scoped Topology

/-!
# §1718 — THE CHEBYSHEV LIMIT FORMS: `ψ(x)/x → 1` and `θ(x)/x → 1`

Completes the classical sextet (ψ/θ/π, each with dlVP error and limit form): dividing §1704/§1705
by `x` gives `|ψ(x)/x − 1| ≤ C·e^{-c√log x} → 0`, same for `θ`.
-/

theorem tendsto_div_one_of_error
    (g : ℝ → ℝ) (C c x₀ : ℝ) (hc : 0 < c) (hx₀ : 1 ≤ x₀)
    (hbound : ∀ x : ℝ, x₀ ≤ x → |g x - x| ≤ C * x * Real.exp (-(c * Real.sqrt (Real.log x)))) :
    Tendsto (fun x : ℝ => g x / x) atTop (𝓝 1) := by
  have hzero : Tendsto (fun x : ℝ => g x / x - 1) atTop (𝓝 0) := by
    apply squeeze_zero_norm' (a := fun x : ℝ => C * Real.exp (-(c * Real.sqrt (Real.log x))))
    · filter_upwards [eventually_ge_atTop (max x₀ 1)] with x hx
      have hxx₀ : x₀ ≤ x := le_trans (le_max_left _ _) hx
      have hx1 : (1 : ℝ) ≤ x := le_trans (le_max_right _ _) hx
      have hx0 : (0 : ℝ) < x := by linarith
      have hb := hbound x hxx₀
      have heq : g x / x - 1 = (g x - x) * (1 / x) := by
        field_simp
      rw [Real.norm_eq_abs, heq, abs_mul, abs_of_pos (by positivity : (0:ℝ) < 1 / x)]
      have hstep : |g x - x| * (1 / x)
          ≤ (C * x * Real.exp (-(c * Real.sqrt (Real.log x)))) * (1 / x) :=
        mul_le_mul_of_nonneg_right hb (by positivity)
      have hcollapse : (C * x * Real.exp (-(c * Real.sqrt (Real.log x)))) * (1 / x)
          = C * Real.exp (-(c * Real.sqrt (Real.log x))) := by
        field_simp
      rw [hcollapse] at hstep
      exact hstep
    · -- C·e^{-c√log x} → 0
      have h1 : Tendsto (fun x : ℝ => c * Real.sqrt (Real.log x)) atTop atTop := by
        apply Tendsto.const_mul_atTop hc
        exact Real.tendsto_sqrt_atTop.comp Real.tendsto_log_atTop
      have h2 : Tendsto (fun x : ℝ => Real.exp (-(c * Real.sqrt (Real.log x)))) atTop (𝓝 0) := by
        have := Real.tendsto_exp_atBot.comp (tendsto_neg_atBot_iff.mpr h1)
        exact this
      have h3 := h2.const_mul C
      simpa using h3
  have hfinal := hzero.add_const 1
  simp only [zero_add] at hfinal
  convert hfinal using 2 with x
  ring

theorem tendsto_psi_div_self : Tendsto (fun x : ℝ => ψ x / x) atTop (𝓝 1) := by
  obtain ⟨C, c, x₀, hC, hc, hx₀, hbound⟩ := pnt_psi_error
  exact tendsto_div_one_of_error ψ C c x₀ hc hx₀ hbound

theorem tendsto_theta_div_self : Tendsto (fun x : ℝ => θ x / x) atTop (𝓝 1) := by
  obtain ⟨C, c, x₀, hC, hc, hx₀, hbound⟩ := pnt_theta_error
  exact tendsto_div_one_of_error θ C c x₀ hc hx₀ hbound

#print axioms tendsto_psi_div_self
#print axioms tendsto_theta_div_self
