import LFunctionOrderNeTop
import LFunctionDivisorFinite
import Mathlib.Analysis.Meromorphic.FactorizedRational
import Mathlib.Analysis.Meromorphic.Divisor

/-!
# §1476 — L divisor factorization on a disk: `LFunction χ = (∏ᶠ divisor)·g`, `g` analytic nonvanishing

The L-analogue of the ζ capstone §1343.  Applying Mathlib's `MeromorphicOn.extract_zeros_poles` to the entire
`LFunction χ` (`χ ≠ 1`) on `ball c R`: it is meromorphic there, has finite order everywhere (§1474, from one
nonvanishing point `z₀`) and finite divisor support (§1475), so it factors as

  `LFunction χ =ᶠ[codiscreteWithin (ball c R)] (∏ᶠ u,(·-u)^{divisor (LFunction χ) u}) • g`

with `g` analytic and nonvanishing on the ball.  This is the divisor factorization that feeds the generic
factor-analyticity bricks (§1383/§1346) and ultimately the Borel–Carathéodory bound §1471 at `F = LFunction χ`,
discharging the `Q` term of the L zero-free region §1470.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Filter

theorem exists_factorization_LFunction {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c : ℂ} {R : ℝ} {z₀ : ℂ} (hz₀ : z₀ ∈ ball c R) (hz₀ne : LFunction χ z₀ ≠ 0) :
    ∃ g : ℂ → ℂ, AnalyticOnNhd ℂ g (ball c R) ∧ (∀ u : ball c R, g u ≠ 0) ∧
      LFunction χ =ᶠ[codiscreteWithin (ball c R)]
        (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R) u)) • g := by
  have hconn : IsConnected (ball c R) := ⟨⟨z₀, hz₀⟩, (convex_ball c R).isPreconnected⟩
  have hana : AnalyticOnNhd ℂ (LFunction χ) Set.univ :=
    (differentiable_LFunction hχ).differentiableOn.analyticOnNhd isOpen_univ
  have hmero : MeromorphicOn (LFunction χ) (ball c R) :=
    hana.meromorphicOn.mono_set (Set.subset_univ _)
  have h₂ : ∀ u : ball c R, meromorphicOrderAt (LFunction χ) u ≠ ⊤ :=
    meromorphicOrderAt_LFunction_ne_top hχ hconn hz₀ hz₀ne
  have h₃ : (MeromorphicOn.divisor (LFunction χ) (ball c R)).support.Finite :=
    divisor_LFunction_ball_support_finite hχ c R
  exact hmero.extract_zeros_poles h₂ h₃

#print axioms exists_factorization_LFunction
