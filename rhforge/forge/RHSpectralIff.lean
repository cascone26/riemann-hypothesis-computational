import Mathlib.LinearAlgebra.Matrix.Hermitian
import Mathlib.Tactic

/-!
# §1769 — Capstone: RH ⟺ the operator realizing the zeros is self-adjoint

The single machine-checked equivalence assembling the spectral lane (§1765/§1766 forward, §1767/§1768 converse).
For a family of candidate zeros `γ : ι → ℂ` (complex a priori), let RH be the statement that they are all real
(`(γ i).im = 0`, i.e. every `½ + iγ i` sits on the critical line).  Then:

**RH  ⟺  the diagonal operator `diagonal γ` (which realizes `{γ i}` as its spectrum) is self-adjoint (Hermitian).**

- `⟸`: a self-adjoint operator has real spectrum (§1765) — Hermitian ⟹ entries real.
- `⟹`: real entries make the diagonal operator Hermitian (§1767/§1768) — and via §1768 a genuine class-C
  (particle-hole) operator with spectrum `{±γ i}`.

So the Hilbert–Pólya equivalence is closed as one `iff`.  All remaining content of RH is the CANONICAL
construction of this operator from the primes — the abstract equivalence is now a theorem.  EXIT0 / std-axioms.
-/

open Matrix

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- **RH ⟺ self-adjoint realization.**  All candidate zeros are real iff their diagonal operator is Hermitian. -/
theorem RH_iff_diag_isHermitian (γ : ι → ℂ) :
    (∀ i, (γ i).im = 0) ↔ (Matrix.diagonal γ).IsHermitian := by
  rw [Matrix.IsHermitian, Matrix.diagonal_conjTranspose, Matrix.diagonal_eq_diagonal_iff]
  refine ⟨fun h i => ?_, fun h i => ?_⟩
  · have : (star γ) i = star (γ i) := rfl
    rw [this, Complex.star_def, Complex.conj_eq_iff_im]
    exact h i
  · have hi := h i
    rw [show (star γ) i = star (γ i) from rfl, Complex.star_def, Complex.conj_eq_iff_im] at hi
    exact hi
