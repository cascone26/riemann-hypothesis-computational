import Mathlib.LinearAlgebra.Matrix.DotProduct
import Mathlib.LinearAlgebra.Matrix.ConjTranspose
import Mathlib.Data.Complex.Basic
import Mathlib.Tactic

set_option maxHeartbeats 1000000

noncomputable section

open Matrix
open scoped ComplexOrder

/-!
# §1728 ★★★★★ — OPERATOR FORM: `Mᴴ M = q·I` ⟹ every eigenvalue has modulus √q

The dimension-free, coordinate-free heart of the critical-circle bound (§1724 was its 2×2
characteristic-polynomial shadow). In Weil's proof the Frobenius endomorphism `F`, under the Rosati
positive involution, satisfies `F† F = q·1` — it is `√q` times an isometry of the (positive-definite)
Rosati inner product. THIS brick is the clean linear-algebra fact behind "eigenvalues on the critical
circle": for ANY complex matrix `M` with `Mᴴ M = q·1` (`q > 0` real), every eigenvalue `λ` satisfies
`normSq λ = q`, i.e. `|λ| = √q`.

Proof: for an eigenvector `v ≠ 0`, `⟨Mv, Mv⟩ = ⟨v, Mᴴ M v⟩ = q⟨v,v⟩` and also `= |λ|²⟨v,v⟩`; since
`⟨v,v⟩ = Σ|vᵢ|² > 0`, cancel to get `|λ|² = q`. The positive-definiteness of the inner product is the
Rosati positivity; the relation `Mᴴ M = q·1` is the arithmetic of Frobenius. Coordinate-free, any
dimension — the true statement of "Frobenius eigenvalues lie on `|λ| = q^{1/2}`" = RH for the curve.
-/

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **Frobenius-as-scaled-isometry ⟹ critical circle.** If `Mᴴ * M = q • 1` with `q > 0`, then every
eigenvalue `λ` of `M` (eigenvector `v ≠ 0`) satisfies `normSq λ = q` — all eigenvalues on `|λ| = √q`. -/
theorem normSq_eigenvalue_of_conjTranspose_mul
    (M : Matrix n n ℂ) (q : ℝ) (hq : 0 < q)
    (hM : Mᴴ * M = (q : ℂ) • 1)
    (l : ℂ) (v : n → ℂ) (hv : v ≠ 0) (heig : M *ᵥ v = l • v) :
    Complex.normSq l = q := by
  have hnv_ne : star v ⬝ᵥ v ≠ 0 := fun h => hv (dotProduct_star_self_eq_zero.mp h)
  -- ⟨Mv, Mv⟩ = q · ⟨v,v⟩
  have key : star (M *ᵥ v) ⬝ᵥ (M *ᵥ v) = (q : ℂ) * (star v ⬝ᵥ v) := by
    rw [Matrix.star_mulVec, ← Matrix.dotProduct_mulVec, Matrix.mulVec_mulVec, hM,
      Matrix.smul_mulVec, Matrix.one_mulVec, dotProduct_smul, smul_eq_mul]
  -- ⟨Mv, Mv⟩ = (star l * l) · ⟨v,v⟩ using Mv = l·v
  have eig : star (M *ᵥ v) ⬝ᵥ (M *ᵥ v) = star l * l * (star v ⬝ᵥ v) := by
    rw [heig, star_smul, dotProduct_smul, dotProduct_comm (star l • star v) v,
      dotProduct_smul, dotProduct_comm v (star v), smul_eq_mul, smul_eq_mul]
    ring
  -- combine and cancel the (nonzero) squared-norm
  have hcancel : star l * l = (q : ℂ) := mul_right_cancel₀ hnv_ne (eig.symm.trans key)
  have hns : (Complex.normSq l : ℂ) = (q : ℂ) := by
    rw [← hcancel, mul_comm, ← starRingEnd_apply, Complex.mul_conj]
  exact_mod_cast hns

/-- Modulus form: `‖λ‖ = √q`. -/
theorem norm_eigenvalue_of_conjTranspose_mul
    (M : Matrix n n ℂ) (q : ℝ) (hq : 0 < q)
    (hM : Mᴴ * M = (q : ℂ) • 1)
    (l : ℂ) (v : n → ℂ) (hv : v ≠ 0) (heig : M *ᵥ v = l • v) :
    ‖l‖ = Real.sqrt q := by
  have hns := normSq_eigenvalue_of_conjTranspose_mul M q hq hM l v hv heig
  have hsq : ‖l‖ ^ 2 = q := by rw [← Complex.normSq_eq_norm_sq]; exact hns
  rw [← hsq, Real.sqrt_sq (norm_nonneg l)]

#print axioms normSq_eigenvalue_of_conjTranspose_mul
#print axioms norm_eigenvalue_of_conjTranspose_mul
