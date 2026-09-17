import Mathlib.LinearAlgebra.Matrix.Hermitian
import Mathlib.Tactic

/-!
# §1775 — Even-symbol Toeplitz matrices are persymmetric (commute with the reversal involution)

The parity ingredient of the Connes–van Suijlekom / Carathéodory-Fejér engine (Connes 2026 program): the
finite Weil quadratic form is a Toeplitz matrix `T i j = c (i - j)` with an EVEN symbol `c (-k) = c k`.  The
reversal `J` (`J i = Fin.rev i = n-1-i`) is an involution, and every such Toeplitz matrix is PERSYMMETRIC:

    `T (rev i) (rev j) = T i j`   (equivalently `J * T * J = T`, i.e. `J` commutes with `T`).

Combined with §1771 (`SimpleEigenParity`: a matrix commuting with an involution has, at a simple eigenvalue,
an eigenvector of definite parity `J v = ± v`), this yields exactly Connes' hypothesis: the minimal eigenvector
of the Toeplitz Weil form is EVEN (palindromic) — the "simple and even" input the CS theorem requires to
conclude the Mellin-transform zeros lie on the critical line.  EXIT0 / std-axioms.
-/

open Matrix

variable {n : ℕ} {𝕜 : Type*} [CommRing 𝕜]

/-- A Toeplitz matrix with symbol `c : ℤ → 𝕜`: entry `(i,j)` depends only on `i - j`. -/
def toeplitz (c : ℤ → 𝕜) : Matrix (Fin n) (Fin n) 𝕜 :=
  fun i j => c ((i : ℤ) - (j : ℤ))

/-- **Persymmetry of an even-symbol Toeplitz matrix.**  If the symbol is even (`c (-k) = c k`), then reversing
both indices leaves the entry unchanged: `T (rev i) (rev j) = T i j`. -/
theorem toeplitz_persymmetric (c : ℤ → 𝕜) (hc : ∀ k : ℤ, c (-k) = c k) (i j : Fin n) :
    toeplitz c (Fin.rev i) (Fin.rev j) = toeplitz c i j := by
  unfold toeplitz
  -- (rev i : ℤ) - (rev j : ℤ) = (n-1-i) - (n-1-j) = j - i = -(i - j)
  have hrev : ((Fin.rev i : Fin n) : ℤ) - ((Fin.rev j : Fin n) : ℤ) = -((i : ℤ) - (j : ℤ)) := by
    simp only [Fin.val_rev]
    have hi : (i : ℕ) < n := i.isLt
    have hj : (j : ℕ) < n := j.isLt
    omega
  rw [hrev, hc]
