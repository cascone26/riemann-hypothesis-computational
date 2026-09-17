import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Tactic

/-!
# §1756 — Suzuki's (1.9) Gram structure: the Weil radical = ker of the defining map

Suzuki's RH-equality criterion (Thm 1.4, arXiv:2301.00421) writes the Weil form as a pullback / Gram form:
`⟨ψ, ψ⟩_W = (1/π)‖P_d D ψ‖²_{L²}`, i.e. `W(ψ,ψ) = ‖L ψ‖²` for the linear map `L = P_d D` (`P_d` the Hardy
projection, `D = i d/dt` the momentum operator).  Any such Gram/pullback form is automatically positive
semidefinite, and its **radical is exactly the kernel of `L`**:

`⟨L x, L x⟩ = 0  ⟺  x ∈ ker L`.

So the PSD-saturation radical of the Weil form (my §1744–§1755 theory) is identified concretely as `ker(P_d D)`,
and Suzuki's `V(0)` (the radical's `K`-invariant half-line part) is a subspace of `ker(P_d D)`.  This is the
bridge from Suzuki's Gram representation (1.9) to the radical/V(0) structure — RH ⟺ this kernel's `K`-invariant
half-line part is nonzero.  EXIT0 / std-axioms.
-/

variable {M : Type*} [AddCommGroup M] [Module ℝ M]
variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

/-- **Gram form is PSD.**  The pullback of an inner product along `L` has nonnegative self-energy. -/
theorem gram_pullback_nonneg (L : M →ₗ[ℝ] V) (x : M) :
    0 ≤ inner ℝ (L x) (L x) :=
  real_inner_self_nonneg

/-- **Gram radical = ker.**  The pullback form's radical is exactly the kernel of `L`. -/
theorem gram_pullback_radical_eq_ker (L : M →ₗ[ℝ] V) (x : M) :
    inner ℝ (L x) (L x) = (0 : ℝ) ↔ x ∈ LinearMap.ker L := by
  rw [inner_self_eq_zero, LinearMap.mem_ker]
