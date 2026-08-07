import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Tactic
/-- **Machine-checked closing inequality of the anchor ∀-theorem (§808).**
Arb-certified constants (m_E ≥ 0.09176, λ₇ ≤ 0.60696, sup|Ω| ≤ 2.0175,
τ_E ≤ 0.0862, B_prime = 0.3579), ε₀ = 1e-4. The worst-case Weil lower bound
Q ≥ m_E(1-s₀) − 2·B_pf·√(s₀(1-s₀)) − C_f·s₀ − ε₀ is POSITIVE (s₀ = ff·ε₀). Uses √(s₀(1-s₀)) ≤ x₀
(rational over-approx, x₀² ≥ s₀(1-s₀)). -/
theorem anchor_eps0_positive :
    (91759/1000000 : ℝ) * (1 - (5/19652 : ℝ)) - 2 * (889314739/400000000 : ℝ) * Real.sqrt ((5/19652 : ℝ) * (1 - (5/19652 : ℝ)))
      - (7912209/5000000 : ℝ) * (5/19652 : ℝ) - (1/10000 : ℝ) ≥ 0 := by
  have hle : Real.sqrt ((5/19652 : ℝ) * (1 - (5/19652 : ℝ))) ≤ (15949/1000000 : ℝ) := by
    rw [show ((15949/1000000 : ℝ) : ℝ) = Real.sqrt ((15949/1000000 : ℝ)^2) by rw [Real.sqrt_sq]; norm_num]
    apply Real.sqrt_le_sqrt; norm_num
  have hnn : Real.sqrt ((5/19652 : ℝ) * (1 - (5/19652 : ℝ))) ≥ 0 := Real.sqrt_nonneg _
  nlinarith [hle, hnn]
