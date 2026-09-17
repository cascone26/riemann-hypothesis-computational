import LogDerivZetaVonMangoldt
import VonMangoldtPolePointwise

/-!
# §1405 — `hpole` for ζ: `Re(-logDeriv ζ σ) ≤ 1/(σ-1) + A` (§1336's `hpole`)

Bridging §1359 (`Re(L Λ σ) ≤ 1/(σ-1) + ‖logDeriv g σ‖`) through §1356 (`-logDeriv ζ = L Λ`) gives the pole bound
in the form §1336's `hpole` needs, with `T₀ = Re(-logDeriv ζ σ)` and `A = ‖logDeriv((s-1)ζ) σ‖`.  EXIT0 / std.
-/

open Complex ArithmeticFunction LSeries DirichletCharacter
open scoped LSeries.notation

theorem hpole_riemannZeta {σ : ℝ} (hσ : 1 < σ) :
    (-logDeriv riemannZeta ((σ : ℂ))).re
      ≤ 1 / (σ - 1) + ‖logDeriv (LFunctionTrivChar₁ 1) ((σ : ℂ))‖ := by
  have h356 : -logDeriv riemannZeta ((σ : ℂ)) = LSeries (↗Λ) ((σ : ℂ)) := by
    rw [logDeriv_riemannZeta_eq_neg_vonMangoldt (by simpa using hσ)]; ring
  rw [h356]
  exact vonMangoldt_pole_bound_pointwise hσ

#print axioms hpole_riemannZeta
