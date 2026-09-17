import DifferentiableOnTwistLimit
import TwistLimitEqLFunction
import LFunctionAnalyticRePos

open Filter Topology DirichletCharacter

/-! §1532 twistLimit χ = LFunction χ on {Re>0} (identity theorem). EXIT0 / std-axioms. -/

theorem twistLimit_eqOn_LFunction {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M} (hχ : χ ≠ 1) :
    Set.EqOn (twistLimit χ) (LFunction χ) {s : ℂ | 0 < s.re} := by
  have hU : IsOpen {s : ℂ | 0 < s.re} := isOpen_lt continuous_const Complex.continuous_re
  have hf : AnalyticOnNhd ℂ (twistLimit χ) {s : ℂ | 0 < s.re} :=
    (differentiableOn_twistLimit hχ).analyticOnNhd hU
  have hg : AnalyticOnNhd ℂ (LFunction χ) {s : ℂ | 0 < s.re} := LFunction_analyticOnNhd_re_pos hχ
  have h2re : (2 : ℂ).re = 2 := by simp
  have hz₀ : (2 : ℂ) ∈ {s : ℂ | 0 < s.re} := by rw [Set.mem_setOf_eq, h2re]; norm_num
  have heq : twistLimit χ =ᶠ[𝓝 (2 : ℂ)] LFunction χ := by
    have hopen1 : {s : ℂ | 1 < s.re} ∈ 𝓝 (2 : ℂ) :=
      (isOpen_lt continuous_const Complex.continuous_re).mem_nhds (by rw [Set.mem_setOf_eq, h2re]; norm_num)
    filter_upwards [hopen1] with s hs
    exact twistLimit_eq_LFunction hχ hs
  exact hf.eqOn_of_preconnected_of_eventuallyEq hg isPreconnected_re_pos hz₀ heq

#print axioms twistLimit_eqOn_LFunction
