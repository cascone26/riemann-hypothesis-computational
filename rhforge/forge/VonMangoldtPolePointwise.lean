import VonMangoldtPole

/-!
# §1359 — pointwise pole bound `Re(L Λ σ) ≤ 1/(σ-1) + ‖logDeriv g σ‖` (§1355's `hpole`)

The pointwise form of §1137, extracted without the `∀ᶠ` wrapper so it is directly consumable at the DLVP test
abscissa `σ = 1 + 1/(2K)`.  From §1136 (`L Λ σ = 1/(σ-1) - logDeriv g σ`, `g = (s-1)ζ` entire, nonzero at 1),
taking real parts: `Re(L Λ σ) = 1/(σ-1) - Re(logDeriv g σ) ≤ 1/(σ-1) + ‖logDeriv g σ‖`.  With §1356
(`Re(L Λ σ) = Re(-logDeriv ζ σ) = T₀`) and `A := ‖logDeriv g σ‖`, this is exactly §1355's `hpole`
(`T₀ ≤ 1/(σ-1) + A`).  EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction LSeries DirichletCharacter
open scoped LSeries.notation

theorem vonMangoldt_pole_bound_pointwise {σ : ℝ} (hσ : 1 < σ) :
    (LSeries (↗Λ) (σ : ℂ)).re
      ≤ 1 / (σ - 1) + ‖logDeriv (LFunctionTrivChar₁ 1) (σ : ℂ)‖ := by
  rw [LSeries_vonMangoldt_eq_pole_sub hσ, Complex.sub_re]
  have hre_pole : (1 / ((σ : ℂ) - 1)).re = 1 / (σ - 1) := by
    rw [show (1 : ℂ) / ((σ : ℂ) - 1) = (((1 / (σ - 1)) : ℝ) : ℂ) by
      push_cast [sub_ne_zero.mpr (show (σ : ℝ) ≠ 1 by linarith)]; ring, Complex.ofReal_re]
  rw [hre_pole]
  have hle : -(logDeriv (LFunctionTrivChar₁ 1) (σ : ℂ)).re
      ≤ ‖logDeriv (LFunctionTrivChar₁ 1) (σ : ℂ)‖ :=
    (neg_le_abs _).trans (Complex.abs_re_le_norm _)
  linarith

#print axioms vonMangoldt_pole_bound_pointwise
