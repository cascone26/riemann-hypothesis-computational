import Mathlib.Analysis.Calculus.LogDeriv
import Mathlib.Analysis.Complex.Basic
import Mathlib.Analysis.Calculus.Deriv.Add

/-!
# §1315 — Log-derivative of the zero-product: `(∏_{ρ∈S}(s-ρ))'/∏(s-ρ) = ∑_{ρ∈S} 1/(s-ρ)`

The "principal part" half of the partial-fractions decomposition of `ζ'/ζ`.  Writing the finite disk-zero
product `P(s) = ∏_{ρ∈S}(s-ρ)`, its logarithmic derivative is the sum of simple poles `∑_{ρ∈S} 1/(s-ρ)`
(Mathlib's `logDeriv_prod`, each factor `s-ρ` contributing `1/(s-ρ)`).  Together with a Borel–Carathéodory
bound on `logDeriv(ζ/P)` (the nonvanishing quotient) this gives the decomposition `ζ'/ζ = ∑1/(s-ρ) + g`
that §1314 turns into `hzero`.  EXIT0 / std-axioms.
-/

open Complex

theorem logDeriv_prod_sub {S : Finset ℂ} {z : ℂ} (hz : ∀ ρ ∈ S, z - ρ ≠ 0) :
    logDeriv (fun s => ∏ ρ ∈ S, (s - ρ)) z = ∑ ρ ∈ S, 1 / (z - ρ) := by
  have hf : ∀ ρ ∈ S, (fun s => s - ρ) z ≠ 0 := fun ρ hρ => hz ρ hρ
  have hd : ∀ ρ ∈ S, DifferentiableAt ℂ (fun s => s - ρ) z :=
    fun ρ _ => (differentiableAt_fun_id).sub_const ρ
  rw [logDeriv_prod hf hd]
  refine Finset.sum_congr rfl (fun ρ _ => ?_)
  have hderiv : deriv (fun s : ℂ => s - ρ) z = 1 := ((hasDerivAt_id z).sub_const ρ).deriv
  rw [logDeriv_apply, hderiv]

#print axioms logDeriv_prod_sub
