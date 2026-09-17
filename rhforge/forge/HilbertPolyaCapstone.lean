import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Tactic

/-!
# §1766 — Hilbert–Pólya capstone: self-adjoint + particle-hole ⟹ RH + functional equation

The spectral-lane capstone, assembling §1763 (PHS pairing), §1765 (self-adjoint ⟹ real):
a single machine-checked theorem stating exactly what the class-C construction would deliver.

Let `H` be self-adjoint (symmetric) on a complex inner-product space, and let `C` implement the particle-hole /
chiral anticommutation `C H = − H C` (`C` injective).  Then for any eigenpair `H v = μ v` (`v ≠ 0`):

1. `μ` is REAL  — the eigenvalue lies on `ℝ`, so the zero `½ + iμ` is on the critical line (RH), and
2. `−μ` is also an eigenvalue, with eigenvector `C v ≠ 0` — the functional equation `ξ(s) = ξ(1−s)` pairing.

So IF the Hilbert–Pólya operator exists as a self-adjoint operator with a particle-hole symmetry and spectrum
`{±γₖ}`, then RH holds and the zeros come in the functional-equation pairs `γ ↔ −γ`.  This is the precise
payoff of the whole class-C construction — every solved RH-analog (Selberg Laplacian, Ihara adjacency, Weil
Frobenius) instantiates this theorem; ζ's open status is exactly the existence of the operator.
EXIT0 / std-axioms.
-/

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℂ E]

open scoped InnerProductSpace

/-- **Hilbert–Pólya capstone.**  Self-adjoint `H` + particle-hole `C H = −H C` ⟹ each eigenvalue is real
(on the critical line) and paired with `−μ` (functional equation), the partner eigenvector `C v` being nonzero. -/
theorem hilbert_polya_rh (H C : E →ₗ[ℂ] E)
    (hsym : ∀ x y, ⟪H x, y⟫_ℂ = ⟪x, H y⟫_ℂ)
    (hchi : C ∘ₗ H = -(H ∘ₗ C)) (hCinj : Function.Injective C)
    (v : E) (hv : v ≠ 0) (μ : ℂ) (heig : H v = μ • v) :
    (starRingEnd ℂ) μ = μ ∧ H (C v) = (-μ) • (C v) ∧ C v ≠ 0 := by
  refine ⟨?_, ?_, ?_⟩
  · -- (1) eigenvalue real (§1765)
    have h1 : ⟪H v, v⟫_ℂ = ⟪v, H v⟫_ℂ := hsym v v
    rw [heig, inner_smul_left, inner_smul_right] at h1
    exact mul_right_cancel₀ (inner_self_ne_zero.mpr hv) h1
  · -- (2) −μ is an eigenvalue with eigenvector C v (§1763)
    have key : μ • C v = -(H (C v)) := by
      have h := LinearMap.congr_fun hchi v
      simp only [LinearMap.comp_apply, LinearMap.neg_apply] at h
      rw [heig, map_smul] at h
      exact h
    have h2 : H (C v) = -(μ • C v) := by rw [key, neg_neg]
    rw [neg_smul]; exact h2
  · -- (3) partner eigenvector is genuine
    intro h; exact hv (hCinj (by rw [h, map_zero]))
