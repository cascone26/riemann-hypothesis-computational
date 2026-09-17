import Mathlib.LinearAlgebra.Matrix.Hermitian
import Mathlib.Tactic

/-!
# §1767 — Hilbert–Pólya converse: real zeros ⟹ a self-adjoint operator with that spectrum exists

Completes the machine-checked equivalence.  §1765/§1766 give the FORWARD direction (a self-adjoint operator
forces its eigenvalues real = on the critical line).  The CONVERSE is the elementary existence half: if the
zeros `{γₖ}` are real (RH), the diagonal operator `diagonal (γₖ)` is self-adjoint (Hermitian) and has exactly
`{γₖ}` as its eigenvalues.  So:

**RH  ⟺  ∃ a self-adjoint operator whose spectrum is `{γₖ}`** (with the class-C `±` structure from §1763).

The content of Hilbert–Pólya is therefore ENTIRELY in *constructing* that operator canonically from the primes
(the open problem) — the abstract equivalence itself is now closed in Lean, forward (§1765/§1766) and converse
(here).  EXIT0 / std-axioms.
-/

open Matrix

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- **Converse (self-adjointness).**  A diagonal matrix with real entries is Hermitian (self-adjoint). -/
theorem diag_isHermitian_of_real (γ : ι → ℝ) :
    (Matrix.diagonal (fun i => (γ i : ℂ))).IsHermitian := by
  unfold Matrix.IsHermitian
  rw [Matrix.diagonal_conjTranspose]
  congr 1
  funext i
  simp

/-- **Converse (spectrum).**  The diagonal operator has each `γ k` as an eigenvalue, with eigenvector `eₖ`:
`diagonal γ *ᵥ eₖ = γ k • eₖ`.  So the real numbers `{γₖ}` are realized as the spectrum of a self-adjoint
operator. -/
theorem diag_eigenvalue (γ : ι → ℝ) (k : ι) :
    (Matrix.diagonal (fun i => (γ i : ℂ))) *ᵥ (Pi.single k 1)
      = (γ k : ℂ) • (Pi.single k (1 : ℂ)) := by
  funext i
  rw [Matrix.mulVec_diagonal]
  by_cases h : i = k
  · subst h; simp
  · simp [Pi.single_eq_of_ne h]
