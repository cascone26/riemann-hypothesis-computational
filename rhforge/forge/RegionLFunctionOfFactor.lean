import LFunctionHdecompHgBOfFactor
import RegionLFunctionConnect
import BetaLtSigmaKLFunction
import Mathlib.NumberTheory.LSeries.Nonvanishing

/-!
# §1500 ★★★ — quantitative L(·,χ) region from the χ factorization data (only the fixed-point remains)

L-analogue of ζ §1425.  Specializes the coordinated decomposition/BC pair §1493 at the DLVP test centre
`c = σK+iγ` (`σK = 1+1/(2K)`, `K = 3A+4B+C`) and threads it through the connector §1498.  Given the §1476
factorization data `(g₀, …)` on `ball c R'`, the sphere bound `Mζ`, the strict factor bound `hgc_lt`, an explicit
`hBbound` (the BC scalar `≤ B`), the pole/growth value bounds, and the kept-zero geometry, any zero `β+iγ` of
`LFunction χ` (χ≠1) obeys `β ≤ 1 - 1/(14K)`.  The `hβ`, `hz₀`, `hz₀ne` are auto-discharged (§1499, self-membership,
nonvanishing on `Re ≥ 1`).  The ONLY thing between this and an unconditional statement is the scalar constant-choice
fixed-point (choosing `A,B,C,R,R',Mζ,γ` self-consistently, `B ≍ log|γ|`), whose arithmetic bricks are already
generic (§1445/§1456/§1458/§1459).  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Filter Classical

theorem region_LFunction_of_factor {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {γ β A B C R R' Mζ : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C) (hK : 0 < 3 * A + 4 * B + C)
    {g₀ : ℂ → ℂ}
    (hga : AnalyticOnNhd ℂ g₀ (ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R'))
    (hgn : ∀ u : ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R', g₀ u ≠ 0)
    (hfac_cod : LFunction χ =ᶠ[codiscreteWithin
        (ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ)
        (ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R') u)) • g₀)
    (hLsphere : ∀ s ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R,
      ‖LFunction χ s‖ ≤ Mζ)
    (hgc_lt : ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ
          (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R').toFinset.filter
        (fun u => ¬ u ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) (R / 2)),
        ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) - u)
          ^ (MeromorphicOn.divisor (LFunction χ)
            (ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R') u).toNat)
        * g₀ (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))‖
      < Mζ / (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ
          (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R').toFinset.filter
        (· ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) (R / 2)),
        (MeromorphicOn.divisor (LFunction χ)
          (ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R') u).toNat))
    (hBbound : 4 * (Real.log (Mζ / (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ
            (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R').toFinset.filter
          (· ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) (R / 2)),
          (MeromorphicOn.divisor (LFunction χ)
            (ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R') u).toNat))
        - Real.log ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ
            (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R').toFinset.filter
          (fun u => ¬ u ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) (R / 2)),
          ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) - u)
            ^ (MeromorphicOn.divisor (LFunction χ)
              (ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R') u).toNat)
          * g₀ (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))‖) / (R / 2) ≤ B)
    (hLρ₀ : LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hρ₀ball : ((β : ℂ) + I * (γ : ℝ))
      ∈ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hp : ((β : ℂ) + I * (γ : ℝ))
      ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) (R / 2))
    (hpole : (-logDeriv riemannZeta (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ))).re
      ≤ 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1) + A)
    (hgrowth : (-logDeriv (LFunction (χ ^ 2))
        (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ))).re ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hcre : 1 < (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)).re := by
    have hpos : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]
    linarith
  obtain ⟨g, hdec, hgBexpr⟩ :=
    LFunction_hdecomp_and_hgB_of_factor hχ hR hRR' hcre hga hgn hfac_cod hLsphere hgc_lt
  have hgB : -B ≤ g.re := le_trans (neg_le_neg hBbound) hgBexpr
  have hβ := beta_lt_sigmaK_of_LFunction_zero hχ hK hLρ₀
  have hz₀ne : LFunction χ (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) ≠ 0 :=
    LFunction_ne_zero_of_one_le_re χ (Or.inl hχ) hcre.le
  exact region_LFunction_connect hχ hK hA hB hC hβ hρ₀ball hLρ₀
    (mem_ball_self (lt_trans hR hRR')) hz₀ne hp hgB hdec hpole hgrowth

#print axioms region_LFunction_of_factor
