import Mathlib.Analysis.Calculus.LogDeriv
import Mathlib.Analysis.Complex.Basic
import Mathlib.Analysis.Calculus.Deriv.ZPow
import Mathlib.Analysis.Calculus.Deriv.Add

/-!
# §1329 — Log-derivative of the divisor product with MULTIPLICITIES: `logDeriv(∏(s-ρ)^{m ρ}) = ∑ m ρ/(z-ρ)`

The multiplicity-aware generalization of §1315 (which used a simple product `∏(s-ρ)`).  With `m : ℂ → ℤ` the
integer multiplicities (matching Mathlib's ℤ-valued `divisor`), `P(s) = ∏_{ρ∈S}(s-ρ)^{m ρ}` has logarithmic
derivative `∑_{ρ∈S} m ρ/(z-ρ)` (Mathlib `logDeriv_fun_zpow`: `logDeriv(f^n) = n·logDeriv f`, each factor giving
`m ρ · 1/(z-ρ)`).  This removes the simple-zero limitation, so the whole region chain works for ζ's zeros with any
multiplicity — the kept zero `ρ₀` then contributes `m ρ₀/(σ-β) ≥ 1/(σ-β)`, only strengthening the bound.
EXIT0 / std-axioms.
-/

open Complex

theorem logDeriv_prod_zpow_sub {S : Finset ℂ} {m : ℂ → ℤ} {z : ℂ} (hz : ∀ ρ ∈ S, z - ρ ≠ 0) :
    logDeriv (fun s => ∏ ρ ∈ S, (s - ρ) ^ (m ρ)) z = ∑ ρ ∈ S, (m ρ : ℂ) / (z - ρ) := by
  rw [logDeriv_prod (f := fun ρ s => (s - ρ) ^ (m ρ))
    (fun ρ hρ => zpow_ne_zero (m ρ) (hz ρ hρ))
    (fun ρ hρ => ((differentiableAt_fun_id).sub_const ρ).zpow (Or.inl (hz ρ hρ)))]
  refine Finset.sum_congr rfl (fun ρ _ => ?_)
  rw [logDeriv_fun_zpow ((differentiableAt_fun_id).sub_const ρ) (m ρ)]
  have hderiv : deriv (fun s : ℂ => s - ρ) z = 1 := ((hasDerivAt_id z).sub_const ρ).deriv
  have hld : logDeriv (fun s : ℂ => s - ρ) z = 1 / (z - ρ) := by
    rw [logDeriv_apply, hderiv]
  rw [hld]; ring

#print axioms logDeriv_prod_zpow_sub
