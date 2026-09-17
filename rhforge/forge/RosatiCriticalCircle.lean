import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Tactic

set_option maxHeartbeats 1000000

noncomputable section

open Complex

/-!
# §1729 ★★★★★ — ROSATI FORM: a √q-isometry of a positive-definite metric ⟹ eigenvalues on |λ|=√q

The basis-free, inner-product-space statement of Weil's RH-for-curves conclusion (§1728 was its matrix
incarnation). In Weil's proof the Rosati involution makes the correspondence pairing a POSITIVE-DEFINITE
Hermitian form `g`, and Frobenius `F` satisfies `g(Fx, Fy) = q · g(x,y)` — it is `√q` times an isometry of
`g`. THIS brick: on any complex inner product space, a linear map `F` with `⟪Fx, Fy⟫ = q⟪x,y⟫` (`q>0`) has
every eigenvalue `λ` on the critical circle `|λ| = √q`.

This is the cleanest form of the template's endpoint: it depends on NOTHING but (i) positive-definiteness
of the metric (Rosati positivity — the arithmetic input) and (ii) the isometry-scaling `F†F = q` (the
Frobenius relation). Both are exactly what a Spec-ℤ Hodge-index object (T-c) must supply; given them, RH
for the associated zeta is immediate. Proof: `q⟪v,v⟫ = ⟪Fv,Fv⟫ = ⟪λv,λv⟫ = |λ|²⟪v,v⟫`, and `⟪v,v⟫ ≠ 0`.
-/

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℂ V]

/-- **Rosati √q-isometry ⟹ critical circle.** If `F` scales the (positive-definite) inner product by
`q > 0` — `⟪F x, F y⟫ = q · ⟪x, y⟫` — then every eigenvalue `λ` of `F` (eigenvector `v ≠ 0`) satisfies
`normSq λ = q`, i.e. `|λ| = √q`. -/
theorem rosati_normSq_eigenvalue
    (F : V →ₗ[ℂ] V) (q : ℝ) (hq : 0 < q)
    (hF : ∀ x y : V, inner ℂ (F x) (F y) = (q : ℂ) * inner ℂ x y)
    (l : ℂ) (v : V) (hv : v ≠ 0) (heig : F v = l • v) :
    Complex.normSq l = q := by
  have hvv : inner ℂ v v ≠ 0 := inner_self_ne_zero.mpr hv
  have e1 : inner ℂ (F v) (F v) = (q : ℂ) * inner ℂ v v := hF v v
  have e2 : inner ℂ (F v) (F v) = (starRingEnd ℂ) l * l * inner ℂ v v := by
    rw [heig, inner_smul_left, inner_smul_right]; ring
  have hcancel : (starRingEnd ℂ) l * l = (q : ℂ) := mul_right_cancel₀ hvv (e2.symm.trans e1)
  have hns : (Complex.normSq l : ℂ) = (q : ℂ) := by
    rw [← Complex.mul_conj l, mul_comm, hcancel]
  exact_mod_cast hns

/-- Modulus form: `‖λ‖ = √q`. -/
theorem rosati_norm_eigenvalue
    (F : V →ₗ[ℂ] V) (q : ℝ) (hq : 0 < q)
    (hF : ∀ x y : V, inner ℂ (F x) (F y) = (q : ℂ) * inner ℂ x y)
    (l : ℂ) (v : V) (hv : v ≠ 0) (heig : F v = l • v) :
    ‖l‖ = Real.sqrt q := by
  have hns := rosati_normSq_eigenvalue F q hq hF l v hv heig
  have hsq : ‖l‖ ^ 2 = q := by rw [← Complex.normSq_eq_norm_sq]; exact hns
  rw [← hsq, Real.sqrt_sq (norm_nonneg l)]

#print axioms rosati_normSq_eigenvalue
#print axioms rosati_norm_eigenvalue
