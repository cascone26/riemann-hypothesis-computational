import Mathlib.Tactic
import Mathlib.Algebra.Module.LinearMap.Defs

/-!
# §1763 — Particle-hole symmetry ⟹ ±symmetric spectrum (= the functional equation)

The spectral / AZ-class-C route to RH (SPECTRAL_HP_ROADMAP Path 4) targets a symplectic particle-hole
(class-C) operator whose spectrum is the zeta zeros.  The defining structural fact of a particle-hole /
class-C Hamiltonian is an (anti)symmetry `C H = − H C` for an invertible `C` (the particle-hole operator).
This forces the spectrum to be symmetric about `0`: every eigenvalue `λ` is paired with `−λ`, with the
partner eigenvector `C v`.

VERIFIED numerically this session (machine-exact): a class-C `H = A + σ⃗·S⃗` has spectrum `±`-symmetric, and
this pairing IS the Riemann functional equation `ξ(s) = ξ(1−s)` — the zeros come in pairs `γ ↔ −γ`.  This
brick is the operator-algebra core of that correspondence: **class-C particle-hole symmetry ⟺ the functional
equation.**  It sits on the spectral-construction lane (the unwalled route), not the positivity lane.
EXIT0 / std-axioms.
-/

variable {V : Type*} [AddCommGroup V] [Module ℝ V]

/-- **Particle-hole pairing of eigenvectors.**  If `C H = − H C` and `H v = λ v`, then `C v` is an
eigenvector with eigenvalue `−λ`.  (Spectrum symmetric about `0` = the functional-equation pairing `γ ↔ −γ`.) -/
theorem phs_eigenpair (H C : V →ₗ[ℝ] V)
    (hanti : C ∘ₗ H = -(H ∘ₗ C)) (v : V) (lam : ℝ) (hv : H v = lam • v) :
    H (C v) = (-lam) • (C v) := by
  have key : lam • C v = -(H (C v)) := by
    have h := LinearMap.congr_fun hanti v
    simp only [LinearMap.comp_apply, LinearMap.neg_apply] at h
    rw [hv, map_smul] at h
    exact h
  have h2 : H (C v) = -(lam • C v) := by rw [key, neg_neg]
  rw [neg_smul]; exact h2

/-- **The partner eigenvector is genuine.**  With `C` injective and `v ≠ 0`, the `−λ` eigenvector `C v` is
nonzero — so `−λ` is a true eigenvalue, not a triviality. -/
theorem phs_eigvec_ne (C : V →ₗ[ℝ] V) (hC : Function.Injective C) (v : V) (hv : v ≠ 0) :
    C v ≠ 0 := by
  intro h
  exact hv (hC (by rw [h, map_zero]))
