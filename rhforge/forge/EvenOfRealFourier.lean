import Mathlib

open Matrix

namespace CMKeyLemma

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **Real Fourier eigenvalue ⟹ even (`γ = F²`).**

If the grading `γ` is the square of the Fourier operator (`γ = F·F`) and `v` is an `F`-eigenvector
with a *real* eigenvalue `c ∈ {1,−1}` (i.e. `F v = ±v`), then `v` is `γ`-even: `γ v = v`.

This is the structural reason CM's approximant `k_λ = ℰ(h_λ)` is even: it is built from `h₀, h₄`, both
`+1` Fourier eigenfunctions (`i⁰ = i⁴ = 1`), hence Fourier-fixed, hence even under the `u ↦ u⁻¹` grading
(`γ = F²`).  Companion to `parity_definite_of_simple`: together they say a *simple, Fourier-real* ground
state is even — reducing CM's even-simple to "the ground state's Fourier eigenvalue is real (±1)".
Machine-verified, axiom-clean. -/
theorem even_of_real_fourier_eigenvector (F γ : Matrix n n ℝ) (v : n → ℝ) (c : ℝ)
    (hγ : γ = F * F) (hc : c = 1 ∨ c = -1) (hFv : F *ᵥ v = c • v) :
    γ *ᵥ v = v := by
  have hstep : γ *ᵥ v = (c * c) • v := by
    rw [hγ, ← Matrix.mulVec_mulVec, hFv, Matrix.mulVec_smul, hFv, smul_smul]
  have hcc : c * c = 1 := by rcases hc with h | h <;> rw [h] <;> ring
  rw [hstep, hcc, one_smul]

#print axioms even_of_real_fourier_eigenvector

end CMKeyLemma
