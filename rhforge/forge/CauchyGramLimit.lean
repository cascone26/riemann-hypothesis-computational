/-
CauchyGramLimit.lean (HP RH center, LEDGER S2919, 2026-11-08)
Closing gaps (a),(b) of S2918:
(b) nonnegative weights w ≥ 0 (zeros with sin(γL/2)=0 allowed): every minor of C·diag(w)·Cᵀ is ≥ 0;
(a) total nonnegativity of minors passes to limits (e.g. truncated zero sums → full zero sum).
-/
import Mathlib
import CauchyGramTP

open Matrix Finset Filter Topology

namespace CauchyGramLim

open CauchyTP

theorem gram_minor_nonneg {n m k : ℕ} (x : Fin n → ℝ) (y : Fin m → ℝ) (w : Fin m → ℝ)
    (hx : StrictMono x) (hy : StrictMono y) (hpos : ∀ i j, 0 < x i + y j)
    (hw : ∀ j, 0 ≤ w j) (r c : Fin k → Fin n) (hr : StrictMono r) (hc : StrictMono c) :
    0 ≤ ((cauchyRect x y * diagonal w * (cauchyRect x y)ᵀ).submatrix r c).det := by
  set C := cauchyRect x y with hC
  have hsplit : (C * diagonal w * Cᵀ).submatrix r c
      = (C.submatrix r id * diagonal w) * (C.submatrix c id)ᵀ := by
    ext i j; simp [Matrix.mul_apply, submatrix_apply, transpose_apply]
  rw [hsplit, CauchyBinet.cauchy_binet]
  refine Finset.sum_nonneg fun S _ => mul_nonneg ?_ ?_
  · have h1 : (C.submatrix r id * diagonal w).submatrix id (CauchyBinet.emb S)
        = C.submatrix r (CauchyBinet.emb S) * diagonal (w ∘ CauchyBinet.emb S) := by
      ext i j; simp [Matrix.mul_diagonal, submatrix_apply]
    rw [h1, det_mul, det_diagonal]
    exact mul_nonneg (minor_pos x y hx hy hpos r _ hr (CauchyBinet.emb S).strictMono).le
      (prod_nonneg fun i _ => hw _)
  · have h2 : ((C.submatrix c id)ᵀ).submatrix (CauchyBinet.emb S) id
        = (C.submatrix c (CauchyBinet.emb S))ᵀ := by ext; rfl
    rw [h2, det_transpose]
    exact (minor_pos x y hx hy hpos c _ hc (CauchyBinet.emb S).strictMono).le

/-- Minor nonnegativity is closed under entrywise limits. -/
theorem minor_nonneg_of_tendsto {n k : ℕ} {ι : Type*} {l : Filter ι} [l.NeBot]
    (G : ι → Matrix (Fin n) (Fin n) ℝ) (G₀ : Matrix (Fin n) (Fin n) ℝ)
    (hG : Tendsto G l (𝓝 G₀)) (r c : Fin k → Fin n)
    (h : ∀ᶠ t in l, 0 ≤ ((G t).submatrix r c).det) :
    0 ≤ (G₀.submatrix r c).det := by
  have hc : Continuous fun M : Matrix (Fin n) (Fin n) ℝ => (M.submatrix r c).det :=
    (continuous_id.matrix_submatrix r c).matrix_det
  exact ge_of_tendsto ((hc.tendsto G₀).comp hG) h

/-- Truncated weighted Cauchy Grams → limit: the limit Gram is totally nonnegative. -/
theorem gram_limit_TN {n k : ℕ} (x : Fin n → ℝ) (hx : StrictMono x)
    (G₀ : Matrix (Fin n) (Fin n) ℝ)
    (M : ℕ → ℕ) (y : ∀ N, Fin (M N) → ℝ) (w : ∀ N, Fin (M N) → ℝ)
    (hy : ∀ N, StrictMono (y N)) (hpos : ∀ N i j, 0 < x i + y N j) (hw : ∀ N j, 0 ≤ w N j)
    (hG : Tendsto (fun N => cauchyRect x (y N) * diagonal (w N) * (cauchyRect x (y N))ᵀ)
      atTop (𝓝 G₀))
    (r c : Fin k → Fin n) (hr : StrictMono r) (hc : StrictMono c) :
    0 ≤ (G₀.submatrix r c).det :=
  minor_nonneg_of_tendsto _ G₀ hG r c (Eventually.of_forall fun N =>
    gram_minor_nonneg x (y N) (w N) hx (hy N) (hpos N) (hw N) r c hr hc)

end CauchyGramLim

#print axioms CauchyGramLim.gram_minor_nonneg
#print axioms CauchyGramLim.minor_nonneg_of_tendsto
#print axioms CauchyGramLim.gram_limit_TN
