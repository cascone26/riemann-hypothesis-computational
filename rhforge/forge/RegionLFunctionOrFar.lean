import RegionLFunctionExplicitScalar
import BetaLtSigmaKLFunction
import HfarHolds

set_option maxHeartbeats 1000000

/-!
# §1544 — L region with the zero-geometry case-split removed (`hp`/`hρ₀ball` gone)

L-analogue of ζ §1436.  CLOSE (`σK-β ≤ R/2`) → §1543 (`hρ₀ball` from `hp`, `hp` itself); FAR (`σK-β > R/2`)
→ `β < σK - R/2 ≤ 1-1/(14K)` via the scalar `hfar` (§1434 for `8 ≤ 7KR`).  `hβ` (`β < σK`) is derived from the
zero via §1499.  Every remaining hypothesis is scalar/geometry or an explicit (σ,γ) log-bound.  EXIT0 / std-axioms.
-/

open Complex Metric DirichletCharacter

theorem region_LFunction_or_far {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) (hχ2 : χ ^ 2 ≠ 1) {γ β A B C R R' : ℝ}
    (hR : 0 < R) (hR1 : R < 1) (hRR' : R < R')
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C) (hK : 0 < 3 * A + 4 * B + C)
    (hfar : (1 + 1 / (2 * (3 * A + 4 * B + C))) - R / 2 ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)))
    (hLρ₀ : LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hBexp : 4 * Real.log ((2 + (N : ℝ)
        + ((1 + 1 / (2 * (3 * A + 4 * B + C))) + |γ| + R) * (N : ℝ)
          / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - R))
        * (3 + 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1))) / (R / 2) ≤ B)
    (hCexp : 4 * Real.log ((2 + (N : ℝ)
        + ((1 + 1 / (2 * (3 * A + 4 * B + C))) + |2 * γ| + R) * (N : ℝ)
          / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - R))
        * (3 + 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1))) / (R / 2) ≤ C)
    (hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ)) : ℂ))‖ ≤ A) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hβ := beta_lt_sigmaK_of_LFunction_zero hχ hK hLρ₀
  by_cases hp : ((β : ℂ) + I * (γ : ℝ))
      ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) (R / 2)
  · have hρ₀ball : ((β : ℂ) + I * (γ : ℝ))
        ∈ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R' :=
      mem_ball.mpr (lt_of_le_of_lt (mem_closedBall.mp hp) (by linarith))
    exact region_LFunction_explicit_scalar hχ hχ2 hR hR1 hRR' hA hB hC hK hρ₀ball hLρ₀ hp
      hBexp hCexp hA_pole
  · rw [mem_closedBall, not_le] at hp
    have hdist : dist ((β : ℂ) + I * (γ : ℝ))
        (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))
        = |β - (1 + 1 / (2 * (3 * A + 4 * B + C)))| := by
      rw [Complex.dist_eq]
      have hsub : ((β : ℂ) + I * (γ : ℝ)) - (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))
          = (((β - (1 + 1 / (2 * (3 * A + 4 * B + C))) : ℝ)) : ℂ) := by push_cast; ring
      rw [hsub, Complex.norm_real, Real.norm_eq_abs]
    rw [hdist] at hp
    have hβlt : β - (1 + 1 / (2 * (3 * A + 4 * B + C))) < 0 := by linarith
    rw [abs_of_neg hβlt] at hp
    linarith

#print axioms region_LFunction_or_far
