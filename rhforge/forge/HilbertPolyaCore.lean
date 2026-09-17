import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Tactic

/-!
# §1765 — The Hilbert–Pólya core: a self-adjoint operator has real eigenvalues (zeros on the line)

The night's synthesis isolated RH to ONE target: build a self-adjoint operator whose spectrum is the zeta
zeros `{γₖ}`.  This brick banks the payoff logic — the reason such an operator would settle RH — as a
machine-checked statement over a complex inner-product space:

**`H` symmetric (self-adjoint) and `H v = μ v` (`v ≠ 0`)  ⟹  `μ` is real.**

If the Hilbert–Pólya operator exists (self-adjoint, spectrum `{γₖ}`), its eigenvalues are automatically real —
i.e. every zero `½ + iγₖ` has `γₖ ∈ ℝ`, so lies on the critical line: RH.  This is exactly why every SOLVED
RH-analog is solved (Selberg's Laplacian, Ihara's adjacency matrix, Weil's Frobenius are each manifestly
self-adjoint / bounded, handing over real spectra for free), and why ζ — lacking such an operator — is open.

Combined with §1763 (particle-hole ⟹ `±μ` pairing = functional equation) and §1764 (no Kramers ⟹ simple), a
self-adjoint class-C operator with spectrum `{±γₖ}` gives: real (RH) + paired (`ξ(s)=ξ(1−s)`) + simple zeros.
This is the exact target of the class-C construction — the operator-existence question is all that remains.
EXIT0 / std-axioms.
-/

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℂ E]

open scoped InnerProductSpace

/-- **Hilbert–Pólya core.**  An eigenvalue of a self-adjoint (symmetric) operator is real. -/
theorem selfadjoint_eigenvalue_real (H : E →ₗ[ℂ] E)
    (hsym : ∀ x y, ⟪H x, y⟫_ℂ = ⟪x, H y⟫_ℂ)
    (v : E) (hv : v ≠ 0) (μ : ℂ) (heig : H v = μ • v) :
    (starRingEnd ℂ) μ = μ := by
  have h1 : ⟪H v, v⟫_ℂ = ⟪v, H v⟫_ℂ := hsym v v
  rw [heig, inner_smul_left, inner_smul_right] at h1
  have hvv : ⟪v, v⟫_ℂ ≠ 0 := inner_self_ne_zero.mpr hv
  exact mul_right_cancel₀ hvv h1
