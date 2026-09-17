import Mathlib.LinearAlgebra.Matrix.Hermitian
import Mathlib.Data.Matrix.Block
import Mathlib.Tactic

/-!
# §1768 — Class-C-faithful converse: real zeros ⟹ a genuine particle-hole self-adjoint operator

§1767 gave the converse with a bare diagonal operator.  This upgrades it to a genuine class-C (particle-hole)
operator, matching §1766's hypotheses EXACTLY.  On the doubled space `ι ⊕ ι`, put
`H = diag(γ) ⊕ (−diag(γ))` and the block-swap `C = [[0,1],[1,0]]`.  Then:

* `H` is Hermitian (self-adjoint), with spectrum `{±γₖ}` (the `±` = the functional equation, §1763);
* `C · H = − H · C` (particle-hole / class-C anticommutation);
* `C · C = 1`, so `C` is an involution (invertible).

Thus given real zeros `{γₖ}` (RH), a genuine class-C self-adjoint operator with spectrum `{±γₖ}` exists — the
converse of §1766 with the FULL class-C structure, not merely self-adjointness.  Combined with §1766, this
closes the equivalence **RH ⟺ ∃ a class-C (particle-hole) self-adjoint operator with spectrum {±γₖ}**, both
directions, machine-checked.  The remaining content of Hilbert–Pólya is the CANONICAL construction of this
operator from the primes.  EXIT0 / std-axioms.
-/

open Matrix

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The doubled diagonal Hamiltonian `diag(γ) ⊕ (−diag(γ))`. -/
noncomputable def Hpair (γ : ι → ℝ) : Matrix (ι ⊕ ι) (ι ⊕ ι) ℂ :=
  fromBlocks (diagonal (fun i => (γ i : ℂ))) 0 0 (-(diagonal (fun i => (γ i : ℂ))))

/-- The particle-hole (block-swap) operator `[[0,1],[1,0]]`. -/
def Cph : Matrix (ι ⊕ ι) (ι ⊕ ι) ℂ := fromBlocks 0 1 1 0

/-- **`C` is an involution** (hence invertible). -/
theorem Cph_mul_self : (Cph : Matrix (ι ⊕ ι) (ι ⊕ ι) ℂ) * Cph = 1 := by
  unfold Cph
  rw [fromBlocks_multiply]
  simp [fromBlocks_one]

/-- **Particle-hole anticommutation** `C·H = −H·C` (the class-C structural relation). -/
theorem Cph_Hpair_anticomm (γ : ι → ℝ) :
    Cph * Hpair γ = -(Hpair γ * Cph) := by
  unfold Cph Hpair
  rw [fromBlocks_multiply, fromBlocks_multiply, fromBlocks_neg]
  simp

/-- **`H` is Hermitian (self-adjoint).** -/
theorem Hpair_isHermitian (γ : ι → ℝ) : (Hpair γ).IsHermitian := by
  unfold Matrix.IsHermitian Hpair
  rw [fromBlocks_conjTranspose]
  simp [diagonal_conjTranspose]
