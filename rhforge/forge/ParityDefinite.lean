import Mathlib

open Matrix

namespace CMKeyLemma

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **γ-simple ⟹ parity-definite eigenvector.**

If `M` commutes with an involution `γ` (`γ² = 1`, `Mγ = γM`) and `v` is an eigenvector for an
eigenvalue `μ` whose eigenspace is one-dimensional (`hsimple`: every `μ`-eigenvector is a scalar
multiple of `v`), then `v` is a `γ`-eigenvector: `γv = v` or `γv = −v`.

This is the rigorous first half of CM's **even-simple**: for the Weil operator, which commutes with the
`u ↦ u⁻¹` grading, a *simple* ground state is automatically parity-definite — so proving "even" reduces to
excluding the odd case.  (The remaining "even, not odd" and "simple" are CM's open analytic input.)
Machine-verified, axiom-clean. -/
theorem parity_definite_of_simple (M γ : Matrix n n ℝ) (v : n → ℝ) (μ : ℝ)
    (hγ2 : γ * γ = 1) (hcomm : M * γ = γ * M)
    (hv : M *ᵥ v = μ • v) (hv0 : v ≠ 0)
    (hsimple : ∀ w : n → ℝ, M *ᵥ w = μ • w → ∃ c : ℝ, w = c • v) :
    γ *ᵥ v = v ∨ γ *ᵥ v = -v := by
  -- γv is also a μ-eigenvector
  have heig : M *ᵥ (γ *ᵥ v) = μ • (γ *ᵥ v) := by
    rw [Matrix.mulVec_mulVec, hcomm, ← Matrix.mulVec_mulVec, hv, Matrix.mulVec_smul]
  obtain ⟨c, hc⟩ := hsimple (γ *ᵥ v) heig      -- γv = c • v
  -- involution ⟹ c² = 1
  have hinv : γ *ᵥ (γ *ᵥ v) = v := by rw [Matrix.mulVec_mulVec, hγ2, Matrix.one_mulVec]
  have hcc : (c * c) • v = v := by
    have h1 : γ *ᵥ (c • v) = v := by rw [← hc]; exact hinv
    rwa [Matrix.mulVec_smul, hc, smul_smul] at h1
  have hz : (c * c - 1) • v = 0 := by rw [sub_smul, one_smul, hcc, sub_self]
  rcases smul_eq_zero.mp hz with hc0 | hbad
  · have hcc1 : c * c = 1 := by linarith [sub_eq_zero.mp hc0]
    rcases mul_self_eq_one_iff.mp hcc1 with h | h
    · left;  rw [hc, h, one_smul]
    · right; rw [hc, h, neg_one_smul]
  · exact absurd hbad hv0

#print axioms parity_definite_of_simple

end CMKeyLemma
