import DiscriminantTraceBound
import Mathlib.LinearAlgebra.BilinearForm.Basic
import Mathlib.Tactic

/-!
# §1743 — Cauchy–Schwarz for a positive-semidefinite symmetric bilinear form

Reusing the discriminant hinge §1726 (`disc_nonpos_of_quadratic_nonneg`): a symmetric bilinear form that is
nonnegative on the diagonal satisfies Cauchy–Schwarz:

`B` symmetric, `∀ x, 0 ≤ B x x`  ⟹  `(B u v)² ≤ (B u u)·(B v v)`.

Proof: `t ↦ B (u + t·v) (u + t·v) = (B v v)·t² + 2(B u v)·t + (B u u) ≥ 0` for all `t`; when `B v v > 0` the
discriminant hinge §1726 gives `(2 B u v)² ≤ 4(B v v)(B u u)`; when `B v v = 0` the linear nonnegativity forces
`B u v = 0`.  This is the positive-definite/semidefinite companion of the Lorentzian index inequality (§1725),
and applies verbatim to the archimedean metric (PSD, in fact PD by §1738).  EXIT0 / std-axioms.
-/

/-- **PSD Cauchy–Schwarz.**  A diagonal-nonnegative symmetric bilinear form obeys `(B u v)² ≤ B u u · B v v`. -/
theorem psd_cauchy_schwarz {M : Type*} [AddCommGroup M] [Module ℝ M]
    (B : LinearMap.BilinForm ℝ M) (hsymm : ∀ x y, B x y = B y x)
    (hpsd : ∀ x, 0 ≤ B x x) (u v : M) :
    (B u v) ^ 2 ≤ (B u u) * (B v v) := by
  have hexp : ∀ t : ℝ, B (u + t • v) (u + t • v)
      = (B v v) * t ^ 2 + 2 * (B u v) * t + (B u u) := by
    intro t
    simp only [map_add, map_smul, LinearMap.add_apply, LinearMap.smul_apply, smul_eq_mul]
    rw [hsymm v u]; ring
  have hq : ∀ t : ℝ, 0 ≤ (B v v) * t ^ 2 + 2 * (B u v) * t + (B u u) := by
    intro t; rw [← hexp t]; exact hpsd _
  rcases lt_or_eq_of_le (hpsd v) with hlt | heq
  · -- B v v > 0 : discriminant hinge
    have hdisc := disc_nonpos_of_quadratic_nonneg (B v v) (2 * (B u v)) (B u u) hlt hq
    nlinarith [hdisc]
  · -- B v v = 0 : linear nonnegativity forces B u v = 0
    have hbv : B v v = 0 := heq.symm
    have hbuv : B u v = 0 := by
      by_contra hne
      have h2ne : 2 * (B u v) ≠ 0 := mul_ne_zero two_ne_zero hne
      have h1 := hq (-(B u u + 1) / (2 * (B u v)))
      rw [hbv] at h1
      simp only [zero_mul, zero_add] at h1
      have hval : 2 * (B u v) * (-(B u u + 1) / (2 * (B u v))) = -(B u u + 1) := by
        field_simp
      rw [hval] at h1
      linarith
    rw [hbuv, hbv]; simp
