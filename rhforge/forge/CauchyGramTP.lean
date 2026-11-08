/-
CauchyGramTP.lean (HP RH center, LEDGER S2918, 2026-11-08)
Total positivity of weighted Cauchy Gram matrices G = C · diag(w) · Cᵀ, C_ij = 1/(x_i + y_j),
x, y strictly increasing, x_i + y_j > 0, w > 0: every k×k minor of G with k ≤ m is POSITIVE.
This is the linear-algebra core of the S2859 band theorem: in the band ν_N < γ₁ the CCM/Shi
zero-side matrix S_mn = Σ_ρ w_ρ/((ν_m+γ)(ν_n+γ)) is exactly such a Gram matrix (RH-conditional
form: all w_ρ = 4 sin²(γL/2)/L ≥ 0 real). Proof: Cauchy–Binet (S2917) + Cauchy TP (S2915).
-/
import Mathlib
import CauchyBinet
import CauchyTP

open Matrix Finset

namespace CauchyGram

open CauchyTP

theorem gram_minor_pos {n m k : ℕ} (x : Fin n → ℝ) (y : Fin m → ℝ) (w : Fin m → ℝ)
    (hx : StrictMono x) (hy : StrictMono y) (hpos : ∀ i j, 0 < x i + y j)
    (hw : ∀ j, 0 < w j) (r c : Fin k → Fin n) (hr : StrictMono r) (hc : StrictMono c)
    (hk : k ≤ m) :
    0 < ((cauchyRect x y * diagonal w * (cauchyRect x y)ᵀ).submatrix r c).det := by
  set C := cauchyRect x y with hC
  have hsplit : (C * diagonal w * Cᵀ).submatrix r c
      = (C.submatrix r id * diagonal w) * (C.submatrix c id)ᵀ := by
    ext i j; simp [Matrix.mul_apply, submatrix_apply, transpose_apply]
  rw [hsplit, CauchyBinet.cauchy_binet]
  apply Finset.sum_pos
  · intro S _
    apply mul_pos
    · have h1 : (C.submatrix r id * diagonal w).submatrix id (CauchyBinet.emb S)
          = C.submatrix r (CauchyBinet.emb S) * diagonal (w ∘ CauchyBinet.emb S) := by
        ext i j; simp [Matrix.mul_diagonal, submatrix_apply]
      rw [h1, det_mul, det_diagonal]
      exact mul_pos (minor_pos x y hx hy hpos r _ hr (CauchyBinet.emb S).strictMono)
        (prod_pos fun i _ => hw _)
    · have h2 : ((C.submatrix c id)ᵀ).submatrix (CauchyBinet.emb S) id
          = (C.submatrix c (CauchyBinet.emb S))ᵀ := by ext; rfl
      rw [h2, det_transpose]
      exact minor_pos x y hx hy hpos c _ hc (CauchyBinet.emb S).strictMono
  · exact ⟨⟨univ.map (Fin.castLEEmb hk), by simp⟩, mem_univ _⟩

/-- Principal case: G is positive definite-in-minors; in particular every diagonal minor > 0. -/
theorem gram_principal_minor_pos {n m k : ℕ} (x : Fin n → ℝ) (y : Fin m → ℝ) (w : Fin m → ℝ)
    (hx : StrictMono x) (hy : StrictMono y) (hpos : ∀ i j, 0 < x i + y j)
    (hw : ∀ j, 0 < w j) (r : Fin k → Fin n) (hr : StrictMono r) (hk : k ≤ m) :
    0 < ((cauchyRect x y * diagonal w * (cauchyRect x y)ᵀ).submatrix r r).det :=
  gram_minor_pos x y w hx hy hpos hw r r hr hr hk

end CauchyGram

#print axioms CauchyGram.gram_minor_pos
#print axioms CauchyGram.gram_principal_minor_pos
