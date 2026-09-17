import SmoothedPNTdlVP
import Psi1RBridge

set_option maxHeartbeats 1000000

open Complex Chebyshev

/-!
# §1704a — THE REAL-FORM SMOOTHED PNT: `|Ψℝ(x) − x²/2| ≤ C·x²·e^{-c√log x}`

Strips the complex packaging off §1703: `psi1C = ofReal ∘ psi1R` (§1607 bridge), the main term
`2π·(x:ℂ)²/2` is the coercion of `2π·x²/2`, so the complex norm is `2π·|Ψℝ(x) − x²/2|` and the
`2π` cancels into the constant.  This is the form the de-smoothing squeeze consumes.
-/

theorem real_smoothed_pnt :
    ∃ C c x₀ : ℝ, 0 < C ∧ 0 < c ∧ 1 ≤ x₀ ∧ ∀ x : ℝ, x₀ ≤ x →
      |psi1R x - x ^ 2 / 2| ≤ C * x ^ 2 * Real.exp (-(c * Real.sqrt (Real.log x))) := by
  obtain ⟨C, c, x₀, hC, hc, hx₀, hbound⟩ := smoothed_pnt_dlvp
  have hπ : (0 : ℝ) < 2 * Real.pi := by positivity
  refine ⟨C / (2 * Real.pi), c, x₀, by positivity, hc, hx₀, ?_⟩
  intro x hx
  have hx0 : (0 : ℝ) < x := by linarith
  have h := hbound x hx
  -- rewrite the complex expression as a real coercion
  have hexpr : (2 * Real.pi : ℝ) • psi1C x - 2 * (Real.pi : ℂ) * ((x : ℂ) ^ (2 : ℂ) / 2)
      = ((2 * Real.pi * (psi1R x - x ^ 2 / 2) : ℝ) : ℂ) := by
    rw [psi1C_eq_ofReal]
    have hcpow : (x : ℂ) ^ (2 : ℂ) = ((x ^ 2 : ℝ) : ℂ) := by
      rw [show (2 : ℂ) = ((2 : ℕ) : ℂ) by norm_num, Complex.cpow_natCast]
      push_cast
      ring
    rw [hcpow]
    rw [Complex.real_smul]
    push_cast
    ring
  rw [hexpr, Complex.norm_real, Real.norm_eq_abs, abs_mul, abs_of_pos hπ] at h
  -- convert the rpow to npow
  have hrpow : x ^ (2 : ℝ) = x ^ 2 := by
    rw [show (2 : ℝ) = ((2 : ℕ) : ℝ) by norm_num, Real.rpow_natCast]
  rw [hrpow] at h
  rw [div_mul_eq_mul_div, div_mul_eq_mul_div, le_div_iff₀ hπ]
  calc |psi1R x - x ^ 2 / 2| * (2 * Real.pi)
      = 2 * Real.pi * |psi1R x - x ^ 2 / 2| := by ring
    _ ≤ C * x ^ 2 * Real.exp (-(c * Real.sqrt (Real.log x))) := h
    _ = C * x ^ 2 * Real.exp (-(c * Real.sqrt (Real.log x))) := rfl

#print axioms real_smoothed_pnt
